"""Surface--surface intersection: the curve where two surfaces cut each other.

``tensorcad.intersect`` answers the triangle question ("do these two facets touch?")
and ``tensorcad.csg`` answers the implicit question ("which side of the field is this
point on?"). Neither produces the thing a B-rep kernel, a trimmed patch or a plate
seam is actually made of: the *curve* along which two surfaces meet.
``GAP-CAD-002`` in ``docs/oss_inventory_cad/GAP-ANCHORS.md`` records that gap, with
the acceptance criterion "cylinder--cylinder arc length within 3% of the analytic
value, every point on both surfaces to 1e-9 L".

**The method, and why this one.** Two surfaces in space meet in a one-dimensional set,
but the problem collapses to two dimensions the moment one of them is parametric. If
``S(u, v)`` is the parametric surface and ``f(x) = 0`` the implicit one, then

    g(u, v) = f(S(u, v))

vanishes exactly on the intersection, so the curve is the zero level set of a single
scalar field on a *flat rectangular parameter domain* -- and that is a problem
marching squares solves with no 3-D geometry at all: sample ``g`` on a grid, find the
grid edges where its sign changes, bisect along each such edge to locate the crossing,
and chain the crossings cell to cell through the edge keys. Every output vertex is a
parameter pair, so it lies on ``S`` to round-off *by construction*; the only residual
left to report is ``f`` at the point.

The alternative -- stepping along ``t = grad(f1) x grad(f2)`` and projecting back onto
both surfaces -- needs a seed per branch, a step size chosen against a curvature nobody
has measured yet, and it stalls on near-tangential contact. The parameter-space route
needs no seeding (the grid is the seed), gets branch and loop structure for free from
the cell topology, and its failure mode is honest: a feature narrower than one cell is
invisible, which the caller can expose by refining ``samples``.

**Periodic parameters are the trap.** A cylinder's ``u`` and a sphere's ``u`` are
angles: ``0`` and ``2 pi`` are the same meridian, so a loop that closes across the seam
would be torn into an open chain and its length would lose the closing chord. Those
surfaces flag ``periodic = (True, False)``; the grid wraps, the seam edge is bisected in
*unwrapped* parameters (``2 pi + eps`` rather than ``eps``), and the chain closes.
``examples/surface_intersection`` includes the seam-crossing case on purpose.

**Length.** A chord sum is short for a curved arc: the chord error per segment behaves
like ``kappa**2 h**3 / 24``, so the relative error is ``O(h**2)``. ``refine_levels``
inserts a midpoint projected back onto ``g = 0`` by a Newton step in parameter space,
halving ``h`` per level, so each level cuts the error by four. Both the raw chord sum
and the refined length are reported, and the example measures the convergence order
rather than assuming it.

**What it does not do.** Tangential contact (a plane grazing a sphere) has no sign
change and is therefore *not* found: the zero set there is a point, not a curve, and
the result is empty -- stated here rather than hidden. Curves that run off the patch
are truncated at the domain boundary and reported as open. Features thinner than the
sampling are missed.

Implicit surfaces are signed distances in this module, so a residual *is* a length and
a threshold like ``1e-9 L`` means what it says. Every ``*_implicit`` factory returns
``|grad| = 1`` (distance to the infinite surface, no caps).
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable, Sequence

Vec3 = tuple[float, float, float]
Vec2 = tuple[float, float]

DEFAULT_SAMPLES = (96, 24)
DEFAULT_REFINE_LEVELS = 2
DEFAULT_TOLERANCE = 1e-12
_BISECTION_STEPS = 90
_NEWTON_STEPS = 8


class SurfaceIntersectionError(ValueError):
    """Raised for input that cannot describe an intersection problem."""


# --------------------------------------------------------------------------- #
# small vector helpers (the package stays standard-library only)
# --------------------------------------------------------------------------- #
def _sub(a: Vec3, b: Vec3) -> Vec3:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _dot(a: Vec3, b: Vec3) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _norm(a: Vec3) -> float:
    return math.sqrt(_dot(a, a))


def _chord(a: Vec3, b: Vec3) -> float:
    return _norm(_sub(a, b))


def _positive(value: float, what: str) -> float:
    number = float(value)
    if not math.isfinite(number):
        raise SurfaceIntersectionError(f"{what} must be finite, got {value!r}")
    if number <= 0.0:
        raise SurfaceIntersectionError(f"{what} must be positive, got {value!r}")
    return number


# --------------------------------------------------------------------------- #
# parametric surfaces
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class PlanePatch:
    """A parallelogram patch: ``S(u, v) = origin + u * du + v * dv`` on ``[0, 1]**2``."""

    origin: Vec3 = (0.0, 0.0, 0.0)
    du: Vec3 = (1.0, 0.0, 0.0)
    dv: Vec3 = (0.0, 1.0, 0.0)
    periodic: tuple[bool, bool] = (False, False)

    def __post_init__(self) -> None:
        for name in ("origin", "du", "dv"):
            if len(getattr(self, name)) != 3:
                raise SurfaceIntersectionError(f"{name} must be a 3-vector")
        if _norm(self.du) <= 0.0 or _norm(self.dv) <= 0.0:
            raise SurfaceIntersectionError("du and dv must be non-degenerate")

    @property
    def domain(self) -> tuple[Vec2, Vec2]:
        return (0.0, 1.0), (0.0, 1.0)

    def evaluate(self, u: float, v: float) -> Vec3:
        return (
            self.origin[0] + u * self.du[0] + v * self.dv[0],
            self.origin[1] + u * self.du[1] + v * self.dv[1],
            self.origin[2] + u * self.du[2] + v * self.dv[2],
        )


@dataclass(frozen=True)
class CylinderPatch:
    """Lateral surface of a ``z``-axis cylinder: ``u`` is the angle, ``v`` is ``z``.

    ``u`` is periodic (``0`` and ``2 pi`` are the same meridian), which is the case
    that tears naive marching apart where a loop crosses its own seam, so it is flagged
    instead of being left to the caller. The patch is a *surface*, not a solid: it has
    no caps, and ``v`` spans ``[centre_z - height / 2, centre_z + height / 2]``.
    """

    radius: float = 1.0
    height: float = 1.0
    centre: Vec3 = (0.0, 0.0, 0.0)
    periodic: tuple[bool, bool] = (True, False)

    def __post_init__(self) -> None:
        _positive(self.radius, "radius")
        _positive(self.height, "height")

    @property
    def domain(self) -> tuple[Vec2, Vec2]:
        half = 0.5 * float(self.height)
        return (0.0, 2.0 * math.pi), (self.centre[2] - half, self.centre[2] + half)

    def evaluate(self, u: float, v: float) -> Vec3:
        radius = float(self.radius)
        return (
            self.centre[0] + radius * math.cos(u),
            self.centre[1] + radius * math.sin(u),
            v,
        )


@dataclass(frozen=True)
class SpherePatch:
    """Latitude/longitude sphere patch: ``u`` is longitude, ``v`` is latitude.

    ``u`` is periodic and ``v`` runs over ``[-pi/2, pi/2]``, so both poles belong to the
    patch: they are singular *parameter* points but every evaluation there still lands
    on the sphere.
    """

    radius: float = 1.0
    centre: Vec3 = (0.0, 0.0, 0.0)
    periodic: tuple[bool, bool] = (True, False)

    def __post_init__(self) -> None:
        _positive(self.radius, "radius")

    @property
    def domain(self) -> tuple[Vec2, Vec2]:
        return (0.0, 2.0 * math.pi), (-0.5 * math.pi, 0.5 * math.pi)

    def evaluate(self, u: float, v: float) -> Vec3:
        radius = float(self.radius)
        cos_v = math.cos(v)
        return (
            self.centre[0] + radius * cos_v * math.cos(u),
            self.centre[1] + radius * cos_v * math.sin(u),
            self.centre[2] + radius * math.sin(v),
        )


# --------------------------------------------------------------------------- #
# implicit surfaces (unit-gradient signed distances)
# --------------------------------------------------------------------------- #
def sphere_implicit(radius: float, centre: Vec3 = (0.0, 0.0, 0.0)) -> Callable[[Vec3], float]:
    """Signed distance ``|x - centre| - radius``: positive outside the sphere."""

    r = _positive(radius, "radius")
    c = (float(centre[0]), float(centre[1]), float(centre[2]))

    def field(point: Vec3) -> float:
        return _chord(point, c) - r

    return field


def cylinder_implicit(
    radius: float, centre: Vec3 = (0.0, 0.0, 0.0), axis: str = "z"
) -> Callable[[Vec3], float]:
    """Signed distance to the *infinite* cylinder about the given coordinate axis."""

    r = _positive(radius, "radius")
    if axis not in ("x", "y", "z"):
        raise SurfaceIntersectionError(f"axis must be 'x', 'y' or 'z', got {axis!r}")
    c = (float(centre[0]), float(centre[1]), float(centre[2]))
    first, second = {"x": (1, 2), "y": (2, 0), "z": (0, 1)}[axis]

    def field(point: Vec3) -> float:
        a = point[first] - c[first]
        b = point[second] - c[second]
        return math.sqrt(a * a + b * b) - r

    return field


def plane_implicit(point: Vec3, normal: Vec3 = (0.0, 0.0, 1.0)) -> Callable[[Vec3], float]:
    """Signed distance to the plane through ``point`` with normal ``normal``."""

    length = _norm(normal)
    if length <= 0.0:
        raise SurfaceIntersectionError("normal must be non-zero")
    unit = (normal[0] / length, normal[1] / length, normal[2] / length)
    origin = (float(point[0]), float(point[1]), float(point[2]))

    def field(probe: Vec3) -> float:
        return _dot(_sub(probe, origin), unit)

    return field


def _as_field(implicit: object) -> Callable[[Vec3], float]:
    """Accept a point callable, or a batch ``signed_distance`` (the package convention).

    :class:`tensorcad.geometry_field.GeometryField` takes a *sequence* of points and
    returns a sequence of values, so it is adapted one point at a time rather than
    evaluated in bulk: the marching asks for single values at parameter pairs.
    """

    distance = getattr(implicit, "signed_distance", None)
    if callable(distance):

        def batched(point: Vec3) -> float:
            values = distance([(point[0], point[1], point[2])])
            if isinstance(values, (int, float)):
                return float(values)
            return float(values[0])

        return batched
    if callable(implicit):
        return implicit  # type: ignore[return-value]
    raise SurfaceIntersectionError(
        "implicit surface must be callable or expose signed_distance(points)"
    )


# --------------------------------------------------------------------------- #
# results
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class IntersectionCurve:
    """One branch of an intersection curve: parameter pairs and the points they give."""

    points: tuple[Vec3, ...]
    params: tuple[Vec2, ...]
    closed: bool
    length: float
    chord_length: float
    max_residual: float

    def __post_init__(self) -> None:
        if len(self.points) != len(self.params):
            raise SurfaceIntersectionError("points and params must have equal length")

    @property
    def closure_defect(self) -> float:
        """Distance between the two ends -- zero for a branch the march closed up."""

        if len(self.points) < 2:
            return 0.0
        return _chord(self.points[0], self.points[-1])


@dataclass(frozen=True)
class IntersectionResult:
    """Everything the march found, plus what it cost and how far off it landed."""

    curves: tuple[IntersectionCurve, ...]
    samples: tuple[int, int]
    refine_levels: int

    @property
    def total_length(self) -> float:
        return sum(curve.length for curve in self.curves)

    @property
    def total_chord_length(self) -> float:
        return sum(curve.chord_length for curve in self.curves)

    @property
    def max_residual(self) -> float:
        return max((curve.max_residual for curve in self.curves), default=0.0)

    @property
    def closed_curves(self) -> tuple[IntersectionCurve, ...]:
        return tuple(curve for curve in self.curves if curve.closed)

    @property
    def open_curves(self) -> tuple[IntersectionCurve, ...]:
        return tuple(curve for curve in self.curves if not curve.closed)

    @property
    def point_count(self) -> int:
        return sum(len(curve.points) for curve in self.curves)


# --------------------------------------------------------------------------- #
# domain plumbing
# --------------------------------------------------------------------------- #
def _domain_of(surface: object) -> tuple[Vec2, Vec2, tuple[bool, bool]]:
    if not hasattr(surface, "evaluate"):
        raise SurfaceIntersectionError("surface must expose evaluate(u, v)")
    ranges = getattr(surface, "domain", None)
    if ranges is None:
        raise SurfaceIntersectionError(
            "surface must expose a domain ((u0, u1), (v0, v1))"
        )
    if callable(ranges):
        ranges = ranges()
    try:
        (u0, u1), (v0, v1) = ranges
        u_range = (float(u0), float(u1))
        v_range = (float(v0), float(v1))
    except (TypeError, ValueError) as error:
        raise SurfaceIntersectionError(
            "surface.domain must be ((u0, u1), (v0, v1))"
        ) from error
    for low, high in (u_range, v_range):
        if not (math.isfinite(low) and math.isfinite(high)) or high <= low:
            raise SurfaceIntersectionError(f"surface domain is empty: {u_range!r}, {v_range!r}")
    flags = tuple(getattr(surface, "periodic", (False, False)))
    if len(flags) != 2:
        raise SurfaceIntersectionError("surface.periodic must be a pair of booleans")
    return u_range, v_range, (bool(flags[0]), bool(flags[1]))


def _wrap(value: float, span: Vec2) -> float:
    period = span[1] - span[0]
    return span[0] + math.fmod(math.fmod(value - span[0], period) + period, period)


def _evaluate(surface: object, u: float, v: float) -> Vec3:
    point = surface.evaluate(u, v)  # type: ignore[attr-defined]
    return (float(point[0]), float(point[1]), float(point[2]))


# --------------------------------------------------------------------------- #
# the march
# --------------------------------------------------------------------------- #
def surface_surface_intersection(
    surface: object,
    implicit: object,
    *,
    samples: Sequence[int] = DEFAULT_SAMPLES,
    refine_levels: int = DEFAULT_REFINE_LEVELS,
    tolerance: float = DEFAULT_TOLERANCE,
) -> IntersectionResult:
    """March the zero set of ``implicit(surface(u, v))`` over the parameter domain.

    ``surface`` needs ``domain`` (``((u0, u1), (v0, v1))``) and ``evaluate(u, v)``, and
    may declare ``periodic`` (a pair of booleans, default ``(False, False)``).
    :class:`tensorcad.nurbs.NURBSSurface` already satisfies the first two, so a NURBS
    patch goes in directly. ``implicit`` is a callable or anything with a
    ``signed_distance`` method (every :class:`tensorcad.geometry_field.GeometryField`
    has one) and is read as a signed distance, so residuals come back as lengths.

    ``samples`` is the number of grid *cells* in ``u`` and ``v``. Raise it until the
    answer stops moving: a feature narrower than one cell is invisible, and that is the
    one failure this method cannot detect on its own. ``tolerance`` is the length below
    which two points count as the same point when branches are classified.
    """

    if len(samples) != 2:
        raise SurfaceIntersectionError("samples must be a pair of cell counts")
    u_cells, v_cells = int(samples[0]), int(samples[1])
    if u_cells < 2 or v_cells < 2:
        raise SurfaceIntersectionError(
            f"samples must be at least (2, 2), got {(u_cells, v_cells)!r}"
        )
    if not 0 <= refine_levels <= 3:
        raise SurfaceIntersectionError(
            f"refine_levels must be between 0 and 3, got {refine_levels!r}"
        )
    if not math.isfinite(tolerance) or tolerance <= 0.0:
        raise SurfaceIntersectionError(f"tolerance must be positive, got {tolerance!r}")

    field = _as_field(implicit)
    u_range, v_range, (periodic_u, periodic_v) = _domain_of(surface)
    du = (u_range[1] - u_range[0]) / u_cells
    dv = (v_range[1] - v_range[0]) / v_cells
    u_nodes = u_cells if periodic_u else u_cells + 1
    v_nodes = v_cells if periodic_v else v_cells + 1

    def g(u: float, v: float) -> float:
        return float(field(_evaluate(surface, u, v)))

    def u_param(index: int) -> float:
        return u_range[0] + index * du

    def v_param(index: int) -> float:
        return v_range[0] + index * dv

    samples_grid: list[list[Vec3]] = [
        [_evaluate(surface, u_param(i), v_param(j)) for j in range(v_nodes)]
        for i in range(u_nodes)
    ]
    values: list[list[float]] = [[float(field(point)) for point in row] for row in samples_grid]

    # g identically zero means the two surfaces *coincide* on the patch: the intersection
    # is the patch itself, not a curve. Left alone, every grid edge would report a
    # crossing and the result would be a cloud of meaningless branches, so this fails
    # loudly instead. The scale is the patch's own bounding-box diagonal, so the test is
    # a relative one and does not depend on the units the caller works in.
    extent = max(
        max(point[axis] for row in samples_grid for point in row)
        - min(point[axis] for row in samples_grid for point in row)
        for axis in range(3)
    )
    peak = max(abs(value) for row in values for value in row)
    if peak <= 1e-12 * max(extent, 1e-12):
        raise SurfaceIntersectionError(
            "implicit field vanishes on the whole patch: the surfaces coincide there, so "
            "the intersection is a surface rather than a curve"
        )

    def sample(i: int, j: int) -> float:
        return values[i % u_nodes if periodic_u else i][j % v_nodes if periodic_v else j]

    def u_next(i: int) -> int:
        return (i + 1) % u_nodes if periodic_u else i + 1

    def v_next(j: int) -> int:
        return (j + 1) % v_nodes if periodic_v else j + 1

    def positive(value: float) -> bool:
        # Exact zeros are folded onto the positive side. A crossing that lands exactly on
        # a grid node is then still caught by the edge leaving it: the two ends of an
        # edge cannot both be non-negative unless the field does not change sign there.
        return value >= 0.0

    # ---- locate a crossing on a grid edge, in unwrapped parameters ------------ #
    crossings: dict[tuple[str, int, int], tuple[Vec2, Vec3]] = {}

    def locate(kind: str, i: int, j: int) -> tuple[Vec2, Vec3]:
        key = (kind, i, j)
        cached = crossings.get(key)
        if cached is not None:
            return cached
        a = (u_param(i), v_param(j))
        if kind == "u":
            # The last u-edge of a periodic grid wraps onto the first meridian: it must be
            # bisected from ``u1`` towards ``u0 + period``, not from ``u1`` back to ``u0``.
            end = u_param(i) + du
            if periodic_u and i + 1 >= u_nodes:
                end = u_range[0] + (u_range[1] - u_range[0])
            b = (end, v_param(j))
        else:
            b = (u_param(i), v_param(j) + dv)
        g_lo = g(*a)
        g_hi = g(*b)
        lo, hi = a, b
        inside = positive(g_lo)
        for _ in range(_BISECTION_STEPS):
            mid = (0.5 * (lo[0] + hi[0]), 0.5 * (lo[1] + hi[1]))
            g_mid = g(*mid)
            if g_mid == 0.0:
                lo = hi = mid
                break
            if positive(g_mid) == inside:
                lo = mid
            else:
                hi = mid
            if abs(hi[0] - lo[0]) <= 4.0 * math.ulp(max(abs(hi[0]), 1.0)) and abs(
                hi[1] - lo[1]
            ) <= 4.0 * math.ulp(max(abs(hi[1]), 1.0)):
                break
        param = (0.5 * (lo[0] + hi[0]), 0.5 * (lo[1] + hi[1]))
        crossing = (param, _evaluate(surface, param[0], param[1]))
        crossings[key] = crossing
        return crossing

    # ---- marching squares over the g = 0 level set ---------------------------- #
    segments: list[tuple[tuple[str, int, int], tuple[str, int, int]]] = []
    for i in range(u_cells):
        i1 = u_next(i)
        for j in range(v_cells):
            j1 = v_next(j)
            s00 = positive(sample(i, j))
            s10 = positive(sample(i1, j))
            s11 = positive(sample(i1, j1))
            s01 = positive(sample(i, j1))
            hits: list[tuple[str, int, int]] = []
            if s00 != s10:
                hits.append(("u", i, j))
            if s10 != s11:
                hits.append(("v", i1, j))
            if s01 != s11:
                hits.append(("u", i, j1))
            if s00 != s01:
                hits.append(("v", i, j))
            if len(hits) == 2:
                segments.append((hits[0], hits[1]))
            elif len(hits) == 4:
                # Saddle: the asymptotic decider, reduced to the cell-centre sign. When
                # the centre shares the sign of corner 00, that corner is connected to
                # itself, so the segments pair (bottom, left) and (right, top).
                centre = g(0.5 * (u_param(i) + u_param(i) + du),
                           0.5 * (v_param(j) + v_param(j) + dv))
                first, second, third, fourth = hits
                if positive(centre) == s00:
                    segments.append((first, fourth))
                    segments.append((second, third))
                else:
                    segments.append((first, second))
                    segments.append((third, fourth))

    # Every edge that carries a crossing is now located once and memoised, so the
    # chaining walk can look points up by edge key instead of recomputing them.
    for left, right in segments:
        locate(*left)
        locate(*right)

    if not segments:
        return IntersectionResult(
            curves=(), samples=(u_cells, v_cells), refine_levels=refine_levels
        )

    chains = _chains(segments)
    curves: list[IntersectionCurve] = []
    for chain in chains:
        curve = _build_curve(
            surface,
            g,
            chain,
            crossings,
            u_range,
            v_range,
            (periodic_u, periodic_v),
            refine_levels,
            tolerance,
        )
        if curve is not None:
            curves.append(curve)
    return IntersectionResult(
        curves=tuple(curves), samples=(u_cells, v_cells), refine_levels=refine_levels
    )


def _chains(
    segments: Sequence[tuple[tuple[str, int, int], tuple[str, int, int]]]
) -> list[list[tuple[str, int, int]]]:
    """Chain cell-local segments through their shared grid edges.

    Segments in adjacent cells name the *same* edge key at the crossing they share, so
    chaining is a dictionary walk and needs no geometric proximity test (which would be
    ambiguous exactly where two branches cross).
    """

    adjacency: dict[tuple[str, int, int], list[int]] = {}
    for index, (left, right) in enumerate(segments):
        adjacency.setdefault(left, []).append(index)
        adjacency.setdefault(right, []).append(index)
    visited: set[int] = set()

    def walk(start: tuple[str, int, int]) -> list[tuple[str, int, int]]:
        chain = [start]
        current = start
        while True:
            options = [index for index in adjacency.get(current, ()) if index not in visited]
            if not options:
                break
            index = options[0]
            visited.add(index)
            left, right = segments[index]
            current = right if left == current else left
            chain.append(current)
            if current == start:
                break
        return chain

    chains: list[list[tuple[str, int, int]]] = []
    for key, links in adjacency.items():
        if len(links) != 1 or links[0] in visited:
            continue
        chain = walk(key)
        if len(chain) >= 2:
            chains.append(chain)
    for key, links in adjacency.items():
        if all(index in visited for index in links):
            continue
        chain = walk(key)
        if len(chain) >= 2:
            chains.append(chain)
    return chains


def _unwrap(
    param: Vec2,
    previous: Vec2,
    u_range: Vec2,
    v_range: Vec2,
    periodic: tuple[bool, bool],
) -> Vec2:
    u, v = param
    if periodic[0]:
        period = u_range[1] - u_range[0]
        while u - previous[0] > 0.5 * period:
            u -= period
        while previous[0] - u > 0.5 * period:
            u += period
    if periodic[1]:
        period = v_range[1] - v_range[0]
        while v - previous[1] > 0.5 * period:
            v -= period
        while previous[1] - v > 0.5 * period:
            v += period
    return (u, v)


def _project(
    g: Callable[[float, float], float],
    param: Vec2,
    u_range: Vec2,
    v_range: Vec2,
    periodic: tuple[bool, bool],
) -> Vec2 | None:
    """Newton step in parameter space, in the direction of ``grad g``, onto ``g = 0``."""

    u, v = param
    step_u = 1e-7 * (u_range[1] - u_range[0])
    step_v = 1e-7 * (v_range[1] - v_range[0])
    start = abs(g(u, v))
    for _ in range(_NEWTON_STEPS):
        value = g(u, v)
        if value == 0.0:
            return (u, v)
        g_u = (g(u + step_u, v) - g(u - step_u, v)) / (2.0 * step_u)
        g_v = (g(u, v + step_v) - g(u, v - step_v)) / (2.0 * step_v)
        denominator = g_u * g_u + g_v * g_v
        if denominator <= 0.0:
            return None
        u -= value * g_u / denominator
        v -= value * g_v / denominator
        if not periodic[0]:
            u = min(max(u, u_range[0]), u_range[1])
        if not periodic[1]:
            v = min(max(v, v_range[0]), v_range[1])
        if abs(value) <= 1e-15 * max(1.0, start):
            return (u, v)
    return (u, v) if abs(g(u, v)) <= 1e-9 * max(1.0, start) else None


def _arc_length(
    surface: object,
    g: Callable[[float, float], float],
    start: Vec2,
    end: Vec2,
    u_range: Vec2,
    v_range: Vec2,
    periodic: tuple[bool, bool],
    levels: int,
) -> float:
    chord = _chord(_evaluate(surface, *start), _evaluate(surface, *end))
    if levels <= 0:
        return chord
    middle = (0.5 * (start[0] + end[0]), 0.5 * (start[1] + end[1]))
    projected = _project(g, middle, u_range, v_range, periodic)
    if projected is None:
        # No better information than the chord: report the chord rather than invent a
        # curvature. This happens where a branch leaves the patch and the midpoint is
        # clamped onto the boundary.
        return chord
    return _arc_length(surface, g, start, projected, u_range, v_range, periodic, levels - 1) + (
        _arc_length(surface, g, projected, end, u_range, v_range, periodic, levels - 1)
    )


def _build_curve(
    surface: object,
    g: Callable[[float, float], float],
    chain: Sequence[tuple[str, int, int]],
    crossings: dict[tuple[str, int, int], tuple[Vec2, Vec3]],
    u_range: Vec2,
    v_range: Vec2,
    periodic: tuple[bool, bool],
    refine_levels: int,
    tolerance: float,
) -> IntersectionCurve | None:
    closed = chain[0] == chain[-1]
    vertices = [crossings[key] for key in chain]

    # Drop repeated vertices -- a crossing found on two edges of one cell at a shared
    # corner otherwise appears twice -- but never drop the closing vertex of a loop:
    # the closing chord is part of its length.
    kept = [vertices[0]]
    for vertex in vertices[1:-1]:
        if _chord(vertex[1], kept[-1][1]) == 0.0:
            continue
        kept.append(vertex)
    if closed:
        kept.append(vertices[-1])
    elif _chord(vertices[-1][1], kept[-1][1]) != 0.0:
        kept.append(vertices[-1])
    if len(kept) < 2:
        return None
    kept_params = [vertex[0] for vertex in kept]
    kept_points = [vertex[1] for vertex in kept]

    # Parameters are unwrapped along the branch so that a loop closing across the periodic
    # seam is refined as one continuous curve: the closing vertex is the first crossing
    # again, so unwrapping carries it a whole period back and the closing chord stays the
    # one-cell step it truly is.
    unwrapped = [kept_params[0]]
    for param in kept_params[1:]:
        unwrapped.append(_unwrap(param, unwrapped[-1], u_range, v_range, periodic))
    length = 0.0
    chord_length = 0.0
    for start, end, point_a, point_b in zip(
        unwrapped[:-1], unwrapped[1:], kept_points[:-1], kept_points[1:]
    ):
        length += _arc_length(surface, g, start, end, u_range, v_range, periodic, refine_levels)
        chord_length += _chord(point_a, point_b)

    # A loop that closes across a *duplicated* seam (a clamped NURBS patch repeats its
    # boundary row) comes back to its start positionally, not through a shared edge key.
    # Close it on position, and add only the seam gap as a chord: re-refining that last
    # step would average two parameters a whole domain apart and walk round the wrong way.
    if not closed and len(kept_points) >= 3:
        extent = max(
            max(abs(point[axis] - kept_points[0][axis]) for point in kept_points)
            for axis in range(3)
        )
        gap = _chord(kept_points[0], kept_points[-1])
        if gap <= tolerance * max(extent, 1.0):
            closed = True
            length += gap
            chord_length += gap

    residual = max(abs(g(*param)) for param in kept_params)
    reported = tuple(
        (
            _wrap(param[0], u_range) if periodic[0] else param[0],
            _wrap(param[1], v_range) if periodic[1] else param[1],
        )
        for param in unwrapped
    )
    return IntersectionCurve(
        points=tuple(kept_points),
        params=reported,
        closed=closed,
        length=length,
        chord_length=chord_length,
        max_residual=residual,
    )