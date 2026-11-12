"""Two independent routes to the thin-walled shear flow, and the differential between them.

``tensorcad.shear_flow`` (merged: one open tree, or one closed loop) and
``tensorcad.shear_flow_multicell`` (a graph: any number of cells) answer the same physical
question by different algebra -- a first-moment sweep, and a linear system whose unknowns are
one constant per wall. The reduction "a multi-cell body whose internal webs are removed must
reproduce the merged solver" is therefore a *test of both*, and it is only meaningful if the
comparison itself is geometry-aware. The first version of that comparison was not, and it
reported three false failures; ``tests/test_shear_flow_multicell.py`` carried them as
``xfail`` items with the note "not wired up correctly in this test". What was wrong is worth
naming, because each is a trap a caller will hit again:

**1. The arc sense is not shared.** A wall is the same wall whichever end you call ``start``;
the merged module sweeps an open tree from its free edges inwards, so its ``start`` is the
junction side, while the graph keeps the caller's declared polylines. Comparing ``q_at(f)``
against ``q_at(f)`` then compares two points mirrored about the wall centre, and the answer
looks wrong by O(1) even when the physics is identical to round-off. Every comparison here
therefore goes through :func:`canonical_q`, which measures ``q`` from the lexicographically
smaller endpoint of the wall -- a convention both routes can be mapped onto without either
being told about it.

**2. A relative tolerance is meaningless at a zero crossing.** ``q`` is zero at a free edge
and at the flange midpoints of a symmetric closed box, and both solvers land there with
round-off noise (``1e-13`` against ``5e-13``). A relative metric scores that as a 100 % error.
The comparators here are **scale-relative**: the error of the comparison is measured against
the section's own flow scale (:func:`flow_scale`), never against the value being compared.

**3. A spliced wall is not an unspliced wall.** ``MultiCellSection`` splits a flange where an
internal web meets it (it must, or the graph is not connected there), while the merged body
has one 4 m flange. So the two chains are compared as *functions of the physical point*
(:func:`field_at`, :func:`difference_on_points`), not wall by wall.

Everything in this module is test scaffolding except :class:`PerimeterOracle`, which is a
deliberate second implementation of the physics (see its docstring). Only the standard
library is required.
"""
from __future__ import annotations

import itertools
import math
from collections.abc import Callable, Iterable, Sequence

Point = tuple[float, float]

#: Coordinates are snapped to this many decimals before two points are called the same one.
_SNAP = 9


def snap(point: Sequence[float]) -> Point:
    return (round(float(point[0]), _SNAP), round(float(point[1]), _SNAP))


# --------------------------------------------------------------------------------------
# Geometry-aware comparison of two flow reports
# --------------------------------------------------------------------------------------


def edges_of(flow) -> tuple:
    """The wall reports of a flow object from either module (``edges`` / ``edge_flows``)."""
    for name in ("edge_flows", "edges"):
        edges = getattr(flow, name, None)
        if edges is not None:
            return tuple(edges)
    raise AssertionError(f"{type(flow).__name__} carries no wall reports")


def canonical_span(edge) -> tuple[Point, Point]:
    """The wall's two endpoints in a fixed order: the lexicographically smaller first."""
    start, end = snap(edge.start), snap(edge.end)
    return (start, end) if start <= end else (end, start)


def canonical_q(edge, fraction: float) -> float:
    """Signed ``q`` at ``fraction`` of the wall measured from its canonical start.

    ``q`` is positive along the canonical direction (smaller endpoint to larger), which is a
    property of the *wall*, not of the solver that reported it. This is the one mapping that
    makes the two modules' reports comparable: neither solver is asked to agree on a
    convention, both are simply re-expressed in a shared one.
    """
    low, _ = canonical_span(edge)
    forward = snap(edge.start) == low
    value = edge.q_at(fraction if forward else 1.0 - fraction)
    return value if forward else -value


def _stations(stations: int) -> tuple[float, ...]:
    return tuple(index / (stations - 1) for index in range(stations))


def flow_scale(flow, stations: int = 21) -> float:
    """The section's own flow scale: the largest ``|q|`` in the report.

    This is the yardstick every tolerance in the reduction tests is divided by. It is one
    number for the whole body, so a wall that happens to carry zero flow is not compared
    against a tolerance proportional to zero.
    """
    biggest = 0.0
    for edge in edges_of(flow):
        for fraction in _stations(stations):
            biggest = max(biggest, abs(edge.q_at(fraction)))
    return biggest


def span_table(flow) -> dict[tuple[Point, Point], object]:
    """``{canonical span: edge}``. Two walls may not share a span."""
    table: dict[tuple[Point, Point], object] = {}
    for edge in edges_of(flow):
        span = canonical_span(edge)
        if span in table:
            raise AssertionError(f"two walls share the span {span}")
        table[span] = edge
    return table


def shared_spans(mine, reference) -> list[tuple[Point, Point]]:
    """The spans present in both bodies, in a stable order."""
    left = span_table(mine)
    right = span_table(reference)
    return sorted(span for span in left if span in right)


def difference_on_spans(
    mine,
    reference,
    spans: Iterable[tuple[Point, Point]] | None = None,
    stations: int = 21,
) -> float:
    """``max |q_mine - q_reference|`` on walls both bodies have, in flow units."""
    left = span_table(mine)
    right = span_table(reference)
    chosen = sorted(spans) if spans is not None else shared_spans(mine, reference)
    if not chosen:
        raise AssertionError("the two bodies have no wall in common")
    worst = 0.0
    for span in chosen:
        if span not in left or span not in right:
            raise AssertionError(f"span {span} is not in both bodies")
        for fraction in _stations(stations):
            worst = max(
                worst,
                abs(canonical_q(left[span], fraction) - canonical_q(right[span], fraction)),
            )
    return worst


def field_at(flow, point: Sequence[float], *, tolerance: float = 1e-9) -> float:
    """Signed ``q`` at a physical point of the section, in its wall's canonical direction."""
    py, pz = float(point[0]), float(point[1])
    if not hasattr(flow, "edge_flows") and not hasattr(flow, "edges"):
        # An oracle-like object reports q directly rather than wall by wall.
        return float(flow.q_at(point))
    found: list[float] = []
    for edge in edges_of(flow):
        (y0, z0), (y1, z1) = snap(edge.start), snap(edge.end)
        length = math.hypot(y1 - y0, z1 - z0)
        if length == 0.0:  # pragma: no cover - zero-length walls are refused upstream
            continue
        along = ((py - y0) * (y1 - y0) + (pz - z0) * (z1 - z0)) / length
        offset = abs((py - y0) * (z1 - z0) - (pz - z0) * (y1 - y0)) / length
        if offset > tolerance or not -tolerance <= along <= length + tolerance:
            continue
        along = min(max(along, 0.0), length)
        forward = snap(edge.start) == canonical_span(edge)[0]
        # ``canonical_q`` wants the fraction measured from the wall's canonical start.
        fraction = along / length if forward else 1.0 - along / length
        found.append(canonical_q(edge, fraction))
    if not found:
        raise AssertionError(f"no wall passes through {(py, pz)}")
    if len(found) > 1 and max(found) - min(found) > 1e-9 * max(1.0, abs(found[0])):
        raise AssertionError(f"the walls through {(py, pz)} disagree: {found}")
    return found[0]


def points_along(walls, *, count: int = 9) -> list[Point]:
    """Physical points on a polyline chain: ``count`` interior stations per straight span.

    Interior stations only, so a point never lands on a junction -- where a spliced body has
    two walls meeting and an unspliced one has a single wall running through.
    """
    points: list[Point] = []
    for wall in walls:
        vertices = [snap(vertex) for vertex in wall.points]
        for (y0, z0), (y1, z1) in itertools.pairwise(vertices):
            for index in range(count):
                fraction = (index + 0.5) / count
                points.append((y0 + fraction * (y1 - y0), z0 + fraction * (z1 - z0)))
    return points


def difference_on_points(mine, reference, points: Sequence[Point]) -> float:
    """``max |q_mine - q_reference|`` at physical points, in flow units."""
    worst = 0.0
    for point in points:
        worst = max(worst, abs(field_at(mine, point) - field_at(reference, point)))
    return worst


def relative_error(measured: float, reference: float, *, scale: float) -> float:
    """Scale-relative error ``|measured - reference| / scale`` -- never divides by a value.

    A relative metric is the wrong yardstick for a field that has zeros in it (a free edge, a
    symmetric flange midpoint): there the *absolute* agreement is round-off and the relative
    one is 100 %. Dividing by the body's flow scale makes the number mean "this fraction of
    the largest flow anywhere in the section", which is what a reader assumes it means.
    """
    if scale <= 0.0:
        raise AssertionError("the flow scale must be positive")
    return abs(measured - reference) / scale


# --------------------------------------------------------------------------------------
# The independent route: integrate the differential relation, never form a first moment.
# --------------------------------------------------------------------------------------


class PerimeterOracle:
    """A closed single cell solved from the *differential* form of the shear-flow relation.

    The merged module computes ``q`` by sweeping a first moment ``Q(s)`` along the section;
    the graph module solves for one constant per wall and imposes no-twist cell by cell. This
    oracle does neither. It integrates

        ``dq/ds = -(V / I) t(s) (v(s) - v_c)``

    along the perimeter from a cut, then fixes the one free constant by Bredt's no-twist
    condition ``CLOSED_INT q / t ds = 0``. Its global sign is fixed the way the modules fix
    theirs -- the resultant must reproduce the applied shear force -- so no sign convention is
    assumed anywhere.

    Because every wall is straight, ``q`` is a quadratic in the arc length on each wall and
    every integral below is evaluated in closed form, so the oracle is round-off accurate
    rather than discretisation limited -- and it shares no algebra with either solver: no
    first moment is formed anywhere, only the differential relation and the no-twist
    integral. ``v_c`` and ``I`` are computed here from the thin-walled definitions.

    The walls are one straight segment each, in the order the caller declares them, forming a
    closed chain; ``direction=-1`` walks the same body the other way round, which must not
    change the answer (a check on the oracle itself).
    """

    def __init__(
        self,
        walls: Sequence,
        *,
        shear_force: float,
        coordinate: Callable[[Point], float],
        cut_wall: int = 0,
        cut_fraction: float = 0.5,
        direction: int = 1,
    ) -> None:
        if direction not in (-1, 1):
            raise ValueError("direction must be +1 (declared order) or -1 (reversed)")
        if not walls:
            raise ValueError("the oracle needs at least one wall")
        for wall in walls:
            if len(wall.points) != 2:
                raise ValueError("the oracle wants one straight wall per wall (split polylines)")
        self.walls = list(walls)
        self.shear_force = float(shear_force)
        self.coordinate = coordinate
        self.direction = direction
        self.centroid, self.second_moment = self._section_properties()
        if self.second_moment <= 0.0:  # pragma: no cover - a degenerate body
            raise ValueError("the section has no second moment about the requested axis")
        self._build(cut_wall, cut_fraction)

    # -- thin-walled section properties, from their definitions -------------------------
    def _section_properties(self) -> tuple[float, float]:
        area = 0.0
        first = 0.0
        for wall in self.walls:
            (y0, z0), (y1, z1) = wall.points
            length = math.hypot(y1 - y0, z1 - z0)
            v0, v1 = self.coordinate((y0, z0)), self.coordinate((y1, z1))
            area += wall.thickness * length
            first += wall.thickness * length * 0.5 * (v0 + v1)
        if area <= 0.0:  # pragma: no cover - a degenerate body
            raise ValueError("the section has no wall area")
        centroid = first / area
        second = 0.0
        for wall in self.walls:
            (y0, z0), (y1, z1) = wall.points
            length = math.hypot(y1 - y0, z1 - z0)
            v0 = self.coordinate((y0, z0)) - centroid
            v1 = self.coordinate((y1, z1)) - centroid
            # INT_0^L (v0 + (v1 - v0) s / L)^2 ds = L (v0^2 + v0 v1 + v1^2) / 3
            second += wall.thickness * length * (v0 * v0 + v0 * v1 + v1 * v1) / 3.0
        return centroid, second

    # -- the integration: closed form on every wall -------------------------------------
    def _build(self, cut_wall: int, cut_fraction: float) -> None:
        if not 0.0 <= cut_fraction <= 1.0:
            raise ValueError("cut_fraction must lie in [0, 1]")
        if not 0 <= cut_wall < len(self.walls):
            raise ValueError(f"cut_wall {cut_wall} is outside the {len(self.walls)} walls")
        # The chain, in the integrating direction, starting and ending at the cut.
        def declared(index: int) -> tuple[Point, Point]:
            wall = self.walls[index]
            p0, p1 = snap(wall.points[0]), snap(wall.points[1])
            return (p0, p1) if self.direction > 0 else (p1, p0)

        cut_start, cut_end = declared(cut_wall)
        cut_thickness = self.walls[cut_wall].thickness
        cut_point = (
            cut_start[0] + cut_fraction * (cut_end[0] - cut_start[0]),
            cut_start[1] + cut_fraction * (cut_end[1] - cut_start[1]),
        )
        chain: list[tuple[Point, Point, float]] = []
        if cut_fraction < 1.0:
            chain.append((cut_point, cut_end, cut_thickness))
        for step in range(1, len(self.walls)):
            index = (cut_wall + step) % len(self.walls) if self.direction > 0 else (
                cut_wall - step) % len(self.walls)
            start, end = declared(index)
            chain.append((start, end, self.walls[index].thickness))
        if cut_fraction > 0.0:
            chain.append((cut_start, cut_point, cut_thickness))
        self._chain = chain

        scale = self.shear_force / self.second_moment
        # Per wall: start point, unit direction, length, thickness, and q_b = A + B s + C s^2
        # with s measured from that wall's start, q_b(0) being the running integral so far.
        records: list[tuple[Point, float, float, float, float, float, float, float]] = []
        running_integral = 0.0  # INT_0^s t (v - v_c) ds' along the walk so far
        arc = 0.0
        for p0, p1, thickness in chain:
            length = math.hypot(p1[0] - p0[0], p1[1] - p0[1])
            if length <= 0.0:  # pragma: no cover - the cut at a vertex makes an empty half
                continue
            v0 = self.coordinate(p0) - self.centroid
            v1 = self.coordinate(p1) - self.centroid
            slope = (v1 - v0) / length
            a = -scale * running_integral
            b = -scale * thickness * v0
            c = -scale * thickness * slope / 2.0
            records.append((p0, length, thickness, a, b, c,
                            (p1[0] - p0[0]) / length, (p1[1] - p0[1]) / length))
            # INT_0^L (running + t (v0 s + slope s^2/2)) ds
            running_integral += thickness * (v0 * length + slope * length * length / 2.0)
            arc += length
        if not records:  # pragma: no cover - a zero-perimeter body
            raise ValueError("the perimeter has no length")
        self.perimeter = arc

        # No twist: CLOSED_INT q / t ds = 0 fixes the free constant. On a wall,
        # INT (A + B s + C s^2) ds = A L + B L^2 / 2 + C L^3 / 3.
        numerator = 0.0
        denominator = 0.0
        for _, length, thickness, a, b, c, _, _ in records:
            numerator += (a * length + b * length ** 2 / 2.0 + c * length ** 3 / 3.0) / thickness
            denominator += length / thickness
        self.circulation = -numerator / denominator

        # Global sign from the resultant, exactly the way the modules fix theirs.
        self._records = records
        self.global_sign = 1
        self.resultant = self._resultant_vector()
        axial = self.resultant[1] if self._axial_is_z() else self.resultant[0]
        self.global_sign = 1 if axial * self.shear_force >= 0.0 else -1

    def _axial_is_z(self) -> bool:
        """Whether the integrated coordinate is ``z`` (a vertical shear)."""
        return self.coordinate((0.0, 1.0)) != self.coordinate((0.0, 0.0))

    def _flow_local(self, record, s: float) -> float:
        _, _, _, a, b, c, _, _ = record
        return self.global_sign * (a + b * s + c * s * s + self.circulation)

    def _resultant_vector(self) -> tuple[float, float]:
        """``(V_y, V_z) = CLOSED_INT q (dy/ds, dz/ds) ds`` -- the applied-force check."""
        vertical = 0.0
        horizontal = 0.0
        for record in self._records:
            _, length, _, a, b, c, dy, dz = record
            moment = (a * length + b * length ** 2 / 2.0 + c * length ** 3 / 3.0
                      + self.circulation * length)
            vertical += self.global_sign * moment * dz
            horizontal += self.global_sign * moment * dy
        return (horizontal, vertical)

    def q_at(self, point: Sequence[float], *, tolerance: float = 1e-9) -> float:
        """``q`` at a physical point on the perimeter, in the wall's canonical direction.

        Canonical is the direction from the wall's lexicographically smaller endpoint, the
        same one :func:`canonical_q` reports for the two solvers -- so this report can be
        compared with theirs point by point without either side being told about a
        convention.
        """
        py, pz = float(point[0]), float(point[1])
        matches: list[float] = []
        for record in self._records:
            p0, length, _, _, _, _, dy, dz = record
            along = (py - p0[0]) * dy + (pz - p0[1]) * dz
            offset = abs((py - p0[0]) * dz - (pz - p0[1]) * dy)
            if offset > tolerance or not -tolerance <= along <= length + tolerance:
                continue
            fraction = min(max(along / length, 0.0), 1.0)
            sense = 1.0 if snap(p0) <= snap((p0[0] + dy * length, p0[1] + dz * length)) else -1.0
            matches.append(sense * self._flow_local(record, fraction * length))
        if not matches:
            raise AssertionError(f"no wall passes through {(py, pz)}")
        if max(matches) - min(matches) > 1e-9 * max(1.0, abs(matches[0])):
            raise AssertionError(f"the oracle disagrees with itself at {(py, pz)}: {matches}")
        return matches[0]
