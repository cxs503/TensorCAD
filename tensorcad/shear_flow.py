"""Thin-walled shear flow: where a shear force actually goes, and why a channel twists.

:mod:`tensorcad.girder` answers *how much bending* a section can carry, and it is explicit
about the question it does not answer: "not a stress calculation ... the section modulus is
the denominator of ``sigma = M / Z``". This module supplies the other half of the
transverse-strength picture that a geometry kernel can own without inventing a load model:

* where the shear force travels in the plating -- the **shear flow** ``q`` (N/m) and the
  **shear stress** ``tau = q / t`` (Pa);
* whether that travel path induces **twist** (the closed-cell circulation ``q0``);
* and where the shear force must be applied so that it does *not* twist the section -- the
  **shear centre**, which is the point a rule check calls the point of twist.

Theory, stated so a consumer can audit the arithmetic
-----------------------------------------------------
The section is a set of **walls**: polylines in the ``(y, z)`` section plane (``y``
transverse starboard-positive, ``z`` vertical up), each with a constant thickness ``t``.
That is the classical thin-walled idealisation: area and first moments ride on the wall
centreline and ``t`` is small compared with the depth.

For a shear force ``V_z`` (the vertical one hull-girder bending produces) the axial stress
is ``sigma_x = M_y z / I_y``, so axial equilibrium of the part of the section cut at arc
length ``s`` requires a shear flow across the cut

    ``q(s) = (V_z / I_y) * Q(s)``,   ``Q(s) = INT_A' (z - z_c) t ds``

with ``A'`` the wall area on one side of the cut and ``I_y = INT (z - z_c)^2 dA`` the
second moment :mod:`tensorcad.girder` already reports as
``second_moment_horizontal``. Two boundary conditions close the problem:

* **open** section (a tree of walls): ``q = 0`` at every free edge, so the flow is swept
  out from the free edges and is fully determined;
* **closed** single-cell section: there is no free edge, so the flow is determined only up
  to a constant circulation ``q0``, fixed by the requirement that pure bending produces *no
  twist*: ``CLOSED_INT q / t ds = 0`` (Bredt compatibility). That circulation is exactly
  why a closed box is torsionally stiff and an open section is not.

The shear centre is then *read off* the flow instead of assumed: the distributed flow is a
force field, its resultant must be the applied shear force ``V``, and the point where that
resultant acts is the shear centre, so ``y_s = c_y + T_z / V_z`` with ``T_z`` the moment of
the flow about the longitudinal axis through the centroid. For a doubly symmetric section
it must land on the centroid; for a channel it must land outside the web by
``e = 3 b^2 / (6 b + h)``, the classical reason a channel twists when loaded through its
centroid.

What is verified, and against what
----------------------------------
``docs/shear_flow.md`` carries the numbers and ``tests/test_shear_flow.py`` recomputes them
(it never reads an artifact). The anchors are analytic; the first four are independent of
this implementation:

1. **Free-edge and junction conditions** -- ``q = 0`` at every free edge and net flow zero
   at every junction, so the flow is solenoidal on the section graph.
2. **Equilibrium** -- the resultant of the flow equals the applied shear force and its
   transverse component vanishes, for open and closed sections alike.
3. **No-twist condition** of the closed cell -- ``CLOSED_INT q / t ds = 0``.
4. **Channel shear centre** ``e = 3 b^2 / (6 b + h)`` from the web centreline, derived in
   the test rather than tabulated.
5. **Open I-beam** -- the web carries the whole shear force and ``q`` follows the textbook
   ``V Q(z) / I`` at every station, recomputed in the test from its own formula.
6. **Doubly symmetric sections** (I, box) put the shear centre on the centroid, which any
   sign error in the flow would break.

Errors are round-off, not discretisation: every integral is evaluated in closed form on
straight walls (the flow is quadratic in the arc length at worst, so three-point Gauss is
exact). That is why the residuals are ``1e-15``-ish rather than percent.

Sign convention
---------------
``q`` is positive in the direction of increasing arc length as reported by
:class:`EdgeFlow` (``start`` -> ``end``): from the junction towards the free edge in an
open section, and around the counter-clockwise loop in a closed one. The global sign of the
field is fixed by the requirement that the flow's resultant equals the applied shear force
-- an equilibrium identity, not a preference -- and the tests assert it. Passing
``shear_force = -1`` reverses every arrow and changes nothing else.

What this module is **not**
---------------------------
* Not a **multi-cell** solver. Two or more closed cells need one circulation per cell
  coupled through the twist compatibility of the shared walls; only the single cell is
  implemented and anything else is rejected rather than answered wrongly (the two-box
  section of a bulk carrier is *not* this).
* Not a section with a **product of inertia**. For ``I_yz != 0`` vertical and horizontal
  flow are coupled through the principal-axis rotation, and that case is rejected
  explicitly instead of being silently approximated.
* Not **shear lag / effective breadth**. The section is fully effective, the same
  assumption :mod:`tensorcad.girder` makes; a wide deck in a slender ship is not fully
  effective near the ends and those corrections belong to the rule set.
* Not a **strength, buckling or fatigue** check and not a load model: no bending moment, no
  material, no allowable stress. ``tau`` is the elastic thin-walled shear stress for the
  shear force the caller supplies; comparing it to a rule is a policy decision.
* Not **plastic** shear, not **warping** stiffness and not a torsional solution:
  ``CLOSED_INT q / t ds = 0`` is the constraint that removes twist, not a twist answer.

Only the Python standard library is required.
"""
from __future__ import annotations

import math
from collections.abc import Iterable, Sequence
from dataclasses import dataclass

Point = tuple[float, float]

__all__ = [
    "AXIS_Y",
    "AXIS_Z",
    "EdgeFlow",
    "MaxShearStress",
    "ShearCentre",
    "ShearFlow",
    "ShearFlowError",
    "ThinWalledSection",
    "Wall",
    "box_section",
    "channel_section",
    "i_section",
    "max_shear_stress",
    "shear_centre",
    "shear_flow",
    "t_section",
]

#: Shear force along ``z``: the hull-girder vertical shear (bending about ``y``).
AXIS_Z = "z"

#: Shear force along ``y``: horizontal shear (bending about ``z``).
AXIS_Y = "y"

#: Vertices closer than this distance are merged into one node (metres).
DEFAULT_MERGE_TOLERANCE = 1e-9

#: Stations reported per wall when the caller does not ask for a count.
DEFAULT_STATIONS = 5

#: Three-point Gauss-Legendre on ``[-1, 1]``, exact for polynomials up to degree five.
_GAUSS3 = (
    (-math.sqrt(0.6), 5.0 / 9.0),
    (0.0, 8.0 / 9.0),
    (math.sqrt(0.6), 5.0 / 9.0),
)


class ShearFlowError(ValueError):
    """Raised when a section cannot be idealised as a thin-walled shear-flow problem."""


# ---------------------------------------------------------------------------
# Exact integrals on one straight wall
# ---------------------------------------------------------------------------

def _int_linear(a0: float, a1: float, length: float) -> float:
    """``INT_0^L a(s) ds`` for ``a`` linear between ``a0`` and ``a1``."""
    return 0.5 * length * (a0 + a1)


def _int_square(a0: float, a1: float, length: float) -> float:
    """``INT_0^L a(s)^2 ds`` for ``a`` linear between ``a0`` and ``a1``."""
    return length * (a0 * a0 + a0 * a1 + a1 * a1) / 3.0


def _int_product(a0: float, a1: float, b0: float, b1: float, length: float) -> float:
    """``INT_0^L a(s) b(s) ds`` for both linear."""
    return length * (a0 * b0 / 3.0 + (a0 * b1 + a1 * b0) / 6.0 + a1 * b1 / 3.0)


def _polyval(coefficients: Sequence[float], s: float) -> float:
    return coefficients[0] + s * (coefficients[1] + s * coefficients[2])


def _poly_integral(coefficients: Sequence[float], length: float) -> float:
    """Exact ``INT_0^L (c0 + c1 s + c2 s^2) ds``."""
    return (
        coefficients[0] * length
        + coefficients[1] * length * length / 2.0
        + coefficients[2] * length ** 3 / 3.0
    )


# ---------------------------------------------------------------------------
# Section description
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Wall:
    """One straight-walled strip: a centreline polyline plus a constant thickness.

    Polyline vertices are section-graph nodes, so two walls connect when they share a
    vertex coordinate, and a wall is subdivided by adding interior vertices (which is how
    it gets split at a junction). ``label`` is carried into the flow report untouched,
    which is how a caller finds "the deck" again.
    """

    points: tuple[Point, ...]
    thickness: float
    label: str = ""

    def __post_init__(self) -> None:
        points = tuple((float(y), float(z)) for y, z in self.points)
        if len(points) < 2:
            raise ShearFlowError(f"wall {self.label!r} needs at least two points")
        for y, z in points:
            if not (math.isfinite(y) and math.isfinite(z)):
                raise ShearFlowError(f"wall {self.label!r} has a non-finite vertex")
        thickness = float(self.thickness)
        if not math.isfinite(thickness) or thickness <= 0.0:
            raise ShearFlowError(
                f"wall {self.label!r} thickness must be finite and positive, "
                f"got {self.thickness!r}"
            )
        object.__setattr__(self, "points", points)
        object.__setattr__(self, "thickness", thickness)


@dataclass(frozen=True)
class _Edge:
    a: int
    b: int
    thickness: float
    label: str


class ThinWalledSection:
    """A tree (open) or a single loop (closed) of thin walls in the ``(y, z)`` plane.

    Construction merges coincident vertices, checks that the walls form exactly one
    connected open tree or exactly one closed loop, and computes the thin-walled area,
    centroid and second moments. Section properties and shear flow are then computed from
    the same idealisation, so they cannot disagree about what the geometry is.
    """

    def __init__(
        self,
        walls: Iterable[Wall],
        *,
        merge_tolerance: float = DEFAULT_MERGE_TOLERANCE,
        closed: bool | None = None,
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

        edges: list[_Edge] = []
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
                edges.append(_Edge(a=first, b=second, thickness=wall.thickness, label=wall.label))

        self.walls = tuple(wall_list)
        self._nodes = tuple(nodes)
        self._edges = tuple(edges)
        self.merge_tolerance = float(merge_tolerance)

        node_count = len(nodes)
        adjacency: list[list[tuple[int, int]]] = [[] for _ in range(node_count)]
        for index, edge in enumerate(edges):
            adjacency[edge.a].append((edge.b, index))
            adjacency[edge.b].append((edge.a, index))
        self._adjacency = tuple(tuple(row) for row in adjacency)

        if any(len(row) == 0 for row in adjacency):
            raise ShearFlowError("the section has an isolated vertex")

        reachable = [False] * node_count
        reachable[0] = True
        stack = [0]
        while stack:
            current = stack.pop()
            for neighbour, _ in adjacency[current]:
                if not reachable[neighbour]:
                    reachable[neighbour] = True
                    stack.append(neighbour)
        if not all(reachable):
            raise ShearFlowError("the walls do not form one connected section")

        degrees = [len(row) for row in adjacency]
        edge_count = len(edges)
        if edge_count == node_count - 1:
            detected = False
        elif edge_count == node_count and all(degree == 2 for degree in degrees):
            detected = True
        else:
            raise ShearFlowError(
                "the walls must form either an open tree (edges = vertices - 1) or exactly "
                "one closed loop (edges = vertices and every vertex of degree two); got "
                f"{edge_count} edges on {node_count} vertices with degrees "
                f"{sorted(set(degrees))}. Multi-cell sections are not supported."
            )
        if closed is not None and bool(closed) != detected:
            raise ShearFlowError(
                f"closed={closed!r} contradicts the wall graph, which is "
                f"{'a closed loop' if detected else 'an open tree'}"
            )
        self.is_closed = detected

        self._compute_properties()

    # -- geometry helpers ---------------------------------------------------
    def _coordinate(self, node: int, axis: str) -> float:
        point = self._nodes[node]
        return point[1] if axis == AXIS_Z else point[0]

    def _length(self, edge: _Edge) -> float:
        (y0, z0), (y1, z1) = self._nodes[edge.a], self._nodes[edge.b]
        return math.hypot(y1 - y0, z1 - z0)

    def _compute_properties(self) -> None:
        area = 0.0
        sum_y = 0.0
        sum_z = 0.0
        for edge in self._edges:
            (y0, z0), (y1, z1) = self._nodes[edge.a], self._nodes[edge.b]
            strip = self._length(edge) * edge.thickness
            area += strip
            sum_y += strip * 0.5 * (y0 + y1)
            sum_z += strip * 0.5 * (z0 + z1)
        if area <= 0.0:
            raise ShearFlowError("the section has no wall area")
        c_y = sum_y / area
        c_z = sum_z / area

        second_z = 0.0  # INT (z - c_z)^2 dA, the bending-about-y inertia
        second_y = 0.0  # INT (y - c_y)^2 dA, the bending-about-z inertia
        product = 0.0
        for edge in self._edges:
            (y0, z0), (y1, z1) = self._nodes[edge.a], self._nodes[edge.b]
            length = self._length(edge) * edge.thickness
            py0, py1 = y0 - c_y, y1 - c_y
            pz0, pz1 = z0 - c_z, z1 - c_z
            second_z += _int_square(pz0, pz1, length)
            second_y += _int_square(py0, py1, length)
            product += _int_product(py0, py1, pz0, pz1, length)

        self.area = area
        self.centroid_y = c_y
        self.centroid_z = c_z
        self.second_moment_horizontal = second_z
        self.second_moment_vertical = second_y
        self.product_of_inertia = product
        self.perimeter = sum(self._length(edge) for edge in self._edges)

    # -- reporting ----------------------------------------------------------
    @property
    def node_count(self) -> int:
        return len(self._nodes)

    @property
    def edge_count(self) -> int:
        return len(self._edges)

    def to_dict(self) -> dict[str, object]:
        return {
            "is_closed": self.is_closed,
            "vertex_count": self.node_count,
            "edge_count": self.edge_count,
            "area": self.area,
            "centroid": [self.centroid_y, self.centroid_z],
            "second_moment_horizontal": self.second_moment_horizontal,
            "second_moment_vertical": self.second_moment_vertical,
            "product_of_inertia": self.product_of_inertia,
            "perimeter": self.perimeter,
            "walls": [
                {
                    "label": wall.label,
                    "thickness": wall.thickness,
                    "points": [list(point) for point in wall.points],
                }
                for wall in self.walls
            ],
        }


# ---------------------------------------------------------------------------
# Convenience sections (the shapes the analytic anchors are written for)
# ---------------------------------------------------------------------------

def _require_positive(**values: float) -> None:
    for name, value in values.items():
        if not math.isfinite(value) or value <= 0.0:
            raise ShearFlowError(f"{name} must be finite and positive, got {value!r}")


def i_section(
    depth: float,
    flange_width: float,
    web_thickness: float,
    flange_thickness: float | None = None,
    *,
    merge_tolerance: float = DEFAULT_MERGE_TOLERANCE,
) -> ThinWalledSection:
    """Doubly symmetric I-section: two horizontal flanges and one vertical web."""
    flange_thickness = web_thickness if flange_thickness is None else flange_thickness
    _require_positive(
        depth=depth,
        flange_width=flange_width,
        web_thickness=web_thickness,
        flange_thickness=flange_thickness,
    )
    half_d, half_b = 0.5 * depth, 0.5 * flange_width
    return ThinWalledSection(
        [
            Wall(((-half_b, half_d), (0.0, half_d), (half_b, half_d)), flange_thickness, "deck"),
            Wall(
                ((-half_b, -half_d), (0.0, -half_d), (half_b, -half_d)),
                flange_thickness,
                "keel",
            ),
            Wall(((0.0, -half_d), (0.0, half_d)), web_thickness, "web"),
        ],
        merge_tolerance=merge_tolerance,
    )


def box_section(
    width: float,
    depth: float,
    web_thickness: float,
    flange_thickness: float | None = None,
    *,
    merge_tolerance: float = DEFAULT_MERGE_TOLERANCE,
) -> ThinWalledSection:
    """Single-cell rectangular box: two vertical webs and two horizontal flanges."""
    flange_thickness = web_thickness if flange_thickness is None else flange_thickness
    _require_positive(
        width=width,
        depth=depth,
        web_thickness=web_thickness,
        flange_thickness=flange_thickness,
    )
    half_d, half_b = 0.5 * depth, 0.5 * width
    corners = [(-half_b, -half_d), (half_b, -half_d), (half_b, half_d), (-half_b, half_d)]
    return ThinWalledSection(
        [
            Wall((corners[0], corners[1]), flange_thickness, "bottom"),
            Wall((corners[1], corners[2]), web_thickness, "starboard"),
            Wall((corners[2], corners[3]), flange_thickness, "top"),
            Wall((corners[3], corners[0]), web_thickness, "port"),
        ],
        merge_tolerance=merge_tolerance,
    )


def channel_section(
    depth: float,
    flange_width: float,
    thickness: float,
    *,
    merge_tolerance: float = DEFAULT_MERGE_TOLERANCE,
) -> ThinWalledSection:
    """Channel: vertical web on ``y = 0``, two flanges projecting towards ``+y``."""
    _require_positive(depth=depth, flange_width=flange_width, thickness=thickness)
    half_d = 0.5 * depth
    return ThinWalledSection(
        [
            Wall(((0.0, -half_d), (0.0, half_d)), thickness, "web"),
            Wall(((0.0, half_d), (flange_width, half_d)), thickness, "top"),
            Wall(((0.0, -half_d), (flange_width, -half_d)), thickness, "bottom"),
        ],
        merge_tolerance=merge_tolerance,
    )


def t_section(
    depth: float,
    flange_width: float,
    thickness: float,
    *,
    merge_tolerance: float = DEFAULT_MERGE_TOLERANCE,
) -> ThinWalledSection:
    """Tee: horizontal flange and a vertical web hanging from its middle."""
    _require_positive(depth=depth, flange_width=flange_width, thickness=thickness)
    half_b = 0.5 * flange_width
    return ThinWalledSection(
        [
            Wall(((-half_b, 0.0), (0.0, 0.0), (half_b, 0.0)), thickness, "flange"),
            Wall(((0.0, 0.0), (0.0, -depth)), thickness, "web"),
        ],
        merge_tolerance=merge_tolerance,
    )

# ---------------------------------------------------------------------------
# Flow construction
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class EdgeFlow:
    """The flow on one wall segment, parameterised by arc length from ``start`` to ``end``.

    ``start`` is the junction end and ``end`` the free-edge end for an open section; for a
    closed section both follow the counter-clockwise loop. ``q`` is positive in that
    direction, so ``start``/``end`` are what make a plot's arrows readable.
    """

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
        """Flow at ``fraction`` of the way from ``start`` to ``end``."""
        if not 0.0 <= fraction <= 1.0:
            raise ShearFlowError("fraction must lie in [0, 1]")
        return _polyval(self.flow_coefficients, fraction * self.length)

    def stress_at(self, fraction: float) -> float:
        """Shear stress ``q / t`` at ``fraction`` of the way along the wall."""
        return self.q_at(fraction) / self.thickness

    def peak_stress(self) -> tuple[float, float]:
        """``(stress, fraction)`` at the stationary point of the quadratic flow, or the end."""
        c0, c1, c2 = self.flow_coefficients
        candidates = [(self.stress_at(0.0), 0.0), (self.stress_at(1.0), 1.0)]
        if c2 != 0.0:
            stationary = -c1 / (2.0 * c2)
            if 0.0 < stationary < self.length:
                candidates.append((_polyval((c0, c1, c2), stationary) / self.thickness,
                                   stationary / self.length))
        return max(candidates, key=lambda item: abs(item[0]))

    def to_dict(self) -> dict[str, object]:
        return {
            "label": self.label,
            "thickness": self.thickness,
            "start": list(self.start),
            "end": list(self.end),
            "length": self.length,
            "s": list(self.s_values),
            "q": list(self.q_values),
            "tau": list(self.stress_values),
        }


@dataclass(frozen=True)
class ShearFlow:
    """A solved flow field plus the residuals that make it auditable.

    The residuals are dimensionless force ratios except ``node_continuity_residual``, which
    is a force in newtons, so a consumer can tell an answer from a bug at a glance.
    """

    axis: str
    shear_force: float
    second_moment: float
    circulation: float
    edges: tuple[EdgeFlow, ...]
    resultant: tuple[float, float]
    equilibrium_residual: float
    transverse_residual: float
    node_continuity_residual: float
    free_edge_residual: float | None
    no_twist_residual: float | None

    @property
    def vertical_resultant(self) -> float:
        return self.resultant[1]

    @property
    def horizontal_resultant(self) -> float:
        return self.resultant[0]

    def edge(self, label: str) -> tuple[EdgeFlow, ...]:
        """Every segment carrying ``label`` (a wall split at a junction appears twice)."""
        return tuple(edge for edge in self.edges if edge.label == label)

    def to_dict(self) -> dict[str, object]:
        return {
            "axis": self.axis,
            "shear_force": self.shear_force,
            "second_moment": self.second_moment,
            "circulation": self.circulation,
            "resultant": list(self.resultant),
            "equilibrium_residual": self.equilibrium_residual,
            "transverse_residual": self.transverse_residual,
            "node_continuity_residual": self.node_continuity_residual,
            "free_edge_residual": self.free_edge_residual,
            "no_twist_residual": self.no_twist_residual,
            "edges": [edge.to_dict() for edge in self.edges],
        }


@dataclass(frozen=True)
class MaxShearStress:
    """Where the section is sheared hardest, and how hard."""

    value: float
    label: str
    fraction: float
    point: Point

    def to_dict(self) -> dict[str, object]:
        return {
            "value": self.value,
            "label": self.label,
            "fraction": self.fraction,
            "point": list(self.point),
        }


@dataclass(frozen=True)
class ShearCentre:
    """The point the shear force must pass through to leave the section untwisted."""

    y: float
    z: float
    torque_vertical: float
    torque_horizontal: float
    vertical_flow: ShearFlow
    horizontal_flow: ShearFlow

    def to_dict(self) -> dict[str, object]:
        return {
            "y": self.y,
            "z": self.z,
            "torque_vertical": self.torque_vertical,
            "torque_horizontal": self.torque_horizontal,
            "vertical_flow": self.vertical_flow.to_dict(),
            "horizontal_flow": self.horizontal_flow.to_dict(),
        }


def _require_uncoupled(section: ThinWalledSection) -> None:
    scale = math.sqrt(
        max(section.second_moment_horizontal, 0.0) * max(section.second_moment_vertical, 0.0)
    )
    if scale <= 0.0:
        raise ShearFlowError("the section has no second moment")
    if abs(section.product_of_inertia) > 1e-9 * scale:
        raise ShearFlowError(
            "the section has a non-zero product of inertia (I_yz != 0): vertical and "
            "horizontal shear flow are coupled and must be solved on the principal axes, "
            "which this module does not do"
        )


@dataclass(frozen=True)
class _Plan:
    """One wall segment's flow description before the ``V / I`` factor and sign are applied.

    ``coefficients`` describe ``Q(s)``, the first moment of the wall area beyond the cut,
    as ``c0 + c1 s + c2 s^2`` with ``s`` measured from ``start`` towards ``end``.
    """

    start: int
    end: int
    coefficients: tuple[float, float, float]
    length: float
    thickness: float
    label: str


def _open_plan(section: ThinWalledSection, axis: str) -> list[_Plan]:
    """Sweep an open tree from its free edges: ``Q`` grows from each free edge inwards."""
    nodes = section._nodes
    edges = section._edges
    adjacency = section._adjacency
    coordinate = section._coordinate
    centroid = section.centroid_z if axis == AXIS_Z else section.centroid_y

    degrees = [len(row) for row in adjacency]
    # The root is a junction (highest degree, ties broken by the lower node index), so every
    # edge's "away from the root" direction runs towards the free edges that feed it.
    root = max(range(len(nodes)), key=lambda index: (degrees[index], -index))

    order = [root]
    children: list[list[tuple[int, int]]] = [[] for _ in nodes]
    visited = [False] * len(nodes)
    visited[root] = True
    cursor = 0
    while cursor < len(order):
        current = order[cursor]
        cursor += 1
        for neighbour, edge_index in adjacency[current]:
            if visited[neighbour]:
                continue
            visited[neighbour] = True
            children[current].append((neighbour, edge_index))
            order.append(neighbour)

    edge_moment = [0.0] * len(edges)
    for index, edge in enumerate(edges):
        v0 = coordinate(edge.a, axis)
        v1 = coordinate(edge.b, axis)
        edge_moment[index] = (
            edge.thickness * section._length(edge) * (0.5 * (v0 + v1) - centroid)
        )

    subtree = [0.0] * len(nodes)
    for node in reversed(order):
        total = 0.0
        for child, edge_index in children[node]:
            total += edge_moment[edge_index] + subtree[child]
        subtree[node] = total

    plans: list[_Plan | None] = [None] * len(edges)
    for node in order:
        for child, edge_index in children[node]:
            edge = edges[edge_index]
            if edge.a == node:
                start, end = edge.a, edge.b
            else:
                start, end = edge.b, edge.a
            length = section._length(edge)
            v0 = coordinate(start, axis)
            v1 = coordinate(end, axis)
            slope = (v1 - v0) / length
            base = 0.5 * (v0 + v1) - centroid  # ((v(s) + v1)/2 - centroid) at s = 0
            half_slope = 0.5 * slope
            thickness = edge.thickness
            # tail(s) = t (L - s) (base + half_slope s), quadratic in s.
            plans[edge_index] = _Plan(
                start=start,
                end=end,
                coefficients=(
                    subtree[child] + thickness * length * base,
                    thickness * (length * half_slope - base),
                    -thickness * half_slope,
                ),
                length=length,
                thickness=thickness,
                label=edge.label,
            )
    if any(plan is None for plan in plans):  # pragma: no cover - a tree reaches every edge
        raise ShearFlowError("internal: an edge was left out of the sweep")
    return [plan for plan in plans if plan is not None]


def _cycle_order(section: ThinWalledSection) -> list[tuple[int, int]]:
    """Walk the loop counter-clockwise (positive centreline area) from node 0."""
    adjacency = section._adjacency
    start = 0
    sequence: list[tuple[int, int]] = []
    previous, current = -1, start
    for _ in range(len(section._edges)):
        options = [(neighbour, index) for neighbour, index in adjacency[current] if neighbour != previous]
        if not options:
            raise ShearFlowError("the closed loop is broken")
        neighbour, _ = options[0]
        sequence.append((current, neighbour))
        previous, current = current, neighbour
        if current == start:
            break
    if current != start or len(sequence) != len(section._edges):
        raise ShearFlowError("the closed loop did not return to its start")
    area2 = 0.0
    for first, second in sequence:
        y0, z0 = section._nodes[first]
        y1, z1 = section._nodes[second]
        area2 += y0 * z1 - y1 * z0
    if area2 < 0.0:
        sequence = [(second, first) for first, second in reversed(sequence)]
    return sequence


def _closed_plan(section: ThinWalledSection, axis: str) -> tuple[list[_Plan], float, float]:
    """Accumulate ``Q`` around the loop.

    Returns the plans, ``CLOSED_INT Q/t ds`` and ``CLOSED_INT ds/t``; the first is what the
    no-twist condition needs to fix the circulation.
    """
    coordinate = section._coordinate
    centroid = section.centroid_z if axis == AXIS_Z else section.centroid_y
    sequence = _cycle_order(section)
    edge_of: dict[tuple[int, int], int] = {}
    for index, edge in enumerate(section._edges):
        edge_of[(min(edge.a, edge.b), max(edge.a, edge.b))] = index

    plans: list[_Plan] = []
    cumulative = 0.0
    twist = 0.0
    length_over_t = 0.0
    for first, second in sequence:
        edge = section._edges[edge_of[(min(first, second), max(first, second))]]
        length = section._length(edge)
        v0 = coordinate(first, axis)
        v1 = coordinate(second, axis)
        slope = (v1 - v0) / length
        thickness = edge.thickness
        coefficients = (cumulative, thickness * (v0 - centroid), thickness * 0.5 * slope)
        twist += _poly_integral(coefficients, length) / thickness
        length_over_t += length / thickness
        plans.append(
            _Plan(
                start=first,
                end=second,
                coefficients=coefficients,
                length=length,
                thickness=thickness,
                label=edge.label,
            )
        )
        cumulative = _polyval(coefficients, length)
    return plans, twist, length_over_t


def _resultant(
    section: ThinWalledSection,
    plans: Sequence[_Plan],
    coefficients: Sequence[Sequence[float]],
) -> tuple[float, float]:
    """Total force carried by the flow field (N), per unit length of the member."""
    nodes = section._nodes
    total_y = 0.0
    total_z = 0.0
    for plan, poly in zip(plans, coefficients):
        y0, z0 = nodes[plan.start]
        y1, z1 = nodes[plan.end]
        dy = (y1 - y0) / plan.length  # unit tangent in the flow direction
        dz = (z1 - z0) / plan.length
        integral = _poly_integral(poly, plan.length)
        total_y += integral * dy
        total_z += integral * dz
    return total_y, total_z


def _torque(
    section: ThinWalledSection,
    plans: Sequence[_Plan],
    coefficients: Sequence[Sequence[float]],
) -> float:
    """Moment of the flow about the longitudinal axis through the centroid (N m per m)."""
    nodes = section._nodes
    c_y, c_z = section.centroid_y, section.centroid_z
    total = 0.0
    for plan, poly in zip(plans, coefficients):
        y0, z0 = nodes[plan.start]
        y1, z1 = nodes[plan.end]
        dy = (y1 - y0) / plan.length
        dz = (z1 - z0) / plan.length
        for node, weight in _GAUSS3:
            s = 0.5 * plan.length * (1.0 + node)
            q = _polyval(poly, s)
            y = y0 + dy * s
            z = z0 + dz * s
            # dT/ds = (y - c_y) q dz/ds - (z - c_z) q dy/ds, cubic in s, so Gauss-3 is exact.
            total += 0.5 * plan.length * weight * q * ((y - c_y) * dz - (z - c_z) * dy)
    return total


def _solve(
    section: ThinWalledSection,
    axis: str,
    shear_force: float,
    stations: int,
) -> tuple[ShearFlow, list[_Plan], list[tuple[float, float, float]]]:
    _require_uncoupled(section)
    if not math.isfinite(shear_force):
        raise ShearFlowError("shear_force must be finite")
    if stations < 2:
        raise ShearFlowError("stations must be at least 2 so that both wall ends appear")

    second_moment = (
        section.second_moment_horizontal if axis == AXIS_Z else section.second_moment_vertical
    )
    if second_moment <= 0.0:
        raise ShearFlowError("the section has no second moment about the requested axis")

    if section.is_closed:
        plans, twist, length_over_t = _closed_plan(section, axis)
        # No twist: CLOSED_INT (q_raw + q0)/t ds = 0  ->  q0 = -CLOSED_INT q_raw/t ds / CLOSED_INT ds/t.
        circulation = -(shear_force / second_moment) * twist / length_over_t
    else:
        plans = _open_plan(section, axis)
        circulation = 0.0

    scale = shear_force / second_moment
    coefficients = [
        (
            scale * plan.coefficients[0] + circulation,
            scale * plan.coefficients[1],
            scale * plan.coefficients[2],
        )
        for plan in plans
    ]
    resultant_y, resultant_z = _resultant(section, plans, coefficients)
    measured = resultant_z if axis == AXIS_Z else resultant_y
    # Sign normalisation: equilibrium fixes the flow's direction only up to one global
    # sign, and the physical choice is the one whose resultant reproduces the applied
    # shear force. This is an equilibrium identity, and the tests assert it.
    if measured * shear_force < 0.0:
        coefficients = [(-poly[0], -poly[1], -poly[2]) for poly in coefficients]
        circulation = -circulation
        resultant_y, resultant_z = _resultant(section, plans, coefficients)

    no_twist: float | None = None
    if section.is_closed:
        # CLOSED_INT q/t ds must vanish for pure bending. It is reported as a dimensionless
        # ratio against the scale it would have if the flow were uniform at its peak.
        integral = 0.0
        length_over_t = 0.0
        peak = 0.0
        for plan, poly in zip(plans, coefficients):
            integral += _poly_integral(poly, plan.length) / plan.thickness
            length_over_t += plan.length / plan.thickness
            for s in (0.0, plan.length):
                peak = max(peak, abs(_polyval(poly, s)))
        denominator = peak * length_over_t
        no_twist = integral / denominator if denominator else 0.0

    edges: list[EdgeFlow] = []
    for plan, poly in zip(plans, coefficients):
        s_values = tuple(plan.length * index / (stations - 1) for index in range(stations))
        q_values = tuple(_polyval(poly, s) for s in s_values)
        edges.append(
            EdgeFlow(
                label=plan.label,
                thickness=plan.thickness,
                start=section._nodes[plan.start],
                end=section._nodes[plan.end],
                length=plan.length,
                flow_coefficients=(poly[0], poly[1], poly[2]),
                s_values=s_values,
                q_values=q_values,
                stress_values=tuple(q / plan.thickness for q in q_values),
            )
        )

    # Junction continuity: at a node the flow that arrives must leave again. A flow is
    # positive from the edge's start to its end, so it leaves the start node and enters the
    # end node; the net outflow of every node must vanish.
    net = [0.0] * len(section._nodes)
    for plan, poly in zip(plans, coefficients):
        net[plan.start] += poly[0]
        net[plan.end] -= _polyval(poly, plan.length)
    node_residual = max(abs(value) for value in net)

    free_edge_residual: float | None = None
    if not section.is_closed:
        worst = 0.0
        for node, row in enumerate(section._adjacency):
            if len(row) != 1:
                continue
            _, edge_index = row[0]
            plan, poly = plans[edge_index], coefficients[edge_index]
            at_start = plan.start == node
            value = poly[0] if at_start else _polyval(poly, plan.length)
            worst = max(worst, abs(value))
        free_edge_residual = worst

    magnitude = abs(shear_force) if shear_force else 1.0
    equilibrium = abs(resultant_z - shear_force) / magnitude if axis == AXIS_Z else abs(
        resultant_y - shear_force
    ) / magnitude
    transverse = abs(resultant_y) / magnitude if axis == AXIS_Z else abs(resultant_z) / magnitude

    flow = ShearFlow(
        axis=axis,
        shear_force=shear_force,
        second_moment=second_moment,
        circulation=circulation,
        edges=tuple(edges),
        resultant=(resultant_y, resultant_z),
        equilibrium_residual=equilibrium,
        transverse_residual=transverse,
        node_continuity_residual=node_residual,
        free_edge_residual=free_edge_residual,
        no_twist_residual=no_twist,
    )
    return flow, plans, coefficients


def shear_flow(
    section: ThinWalledSection,
    shear_force: float,
    *,
    axis: str = AXIS_Z,
    stations: int = DEFAULT_STATIONS,
) -> ShearFlow:
    """Solve the thin-walled shear flow for one shear force.

    ``shear_force`` is in newtons and acts along ``axis`` (``"z"`` vertical by default);
    ``stations`` is how many points per wall the report carries (the flow is a quadratic
    everywhere, so five is already exact -- more stations are for plotting, not accuracy).
    """
    if axis not in (AXIS_Y, AXIS_Z):
        raise ShearFlowError(f"axis must be {AXIS_Z!r} or {AXIS_Y!r}, got {axis!r}")
    flow, _, _ = _solve(section, axis, shear_force, stations)
    return flow


def shear_centre(section: ThinWalledSection, *, unit_force: float = 1.0) -> ShearCentre:
    """Locate the shear centre: the point a shear force must pass through to cause no twist.

    Both directions are solved (vertical shear about the horizontal axis, horizontal shear
    about the vertical one) because the shear centre is a point, not a distance. For a
    section that is symmetric about the vertical axis the two solutions must agree with the
    analytic symmetry axis; a section with ``I_yz != 0`` is rejected.
    """
    _require_uncoupled(section)
    if not math.isfinite(unit_force) or unit_force == 0.0:
        raise ShearFlowError("unit_force must be finite and non-zero")

    vertical_flow, vertical_plans, vertical_coefficients = _solve(
        section, AXIS_Z, unit_force, DEFAULT_STATIONS
    )
    horizontal_flow, horizontal_plans, horizontal_coefficients = _solve(
        section, AXIS_Y, unit_force, DEFAULT_STATIONS
    )
    torque_vertical = _torque(section, vertical_plans, vertical_coefficients)
    torque_horizontal = _torque(section, horizontal_plans, horizontal_coefficients)
    return ShearCentre(
        y=section.centroid_y + torque_vertical / unit_force,
        z=section.centroid_z - torque_horizontal / unit_force,
        torque_vertical=torque_vertical,
        torque_horizontal=torque_horizontal,
        vertical_flow=vertical_flow,
        horizontal_flow=horizontal_flow,
    )


def max_shear_stress(
    section: ThinWalledSection,
    shear_force: float,
    *,
    axis: str = AXIS_Z,
) -> MaxShearStress:
    """Largest ``|tau| = |q| / t`` in the section and where it sits."""
    flow = shear_flow(section, shear_force, axis=axis)
    best: MaxShearStress | None = None
    for edge in flow.edges:
        value, fraction = edge.peak_stress()
        y0, z0 = edge.start
        y1, z1 = edge.end
        point = (y0 + fraction * (y1 - y0), z0 + fraction * (z1 - z0))
        candidate = MaxShearStress(
            value=value, label=edge.label, fraction=fraction, point=point
        )
        if best is None or abs(candidate.value) > abs(best.value):
            best = candidate
    assert best is not None
    return best
