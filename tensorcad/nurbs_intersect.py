"""NURBS x NURBS intersection: the curve where two *parametric* patches meet, plus the
trim domain that curve cuts on each of them.

``tensorcad.surface_intersection`` answers one half of the question -- "where does a
parametric surface meet an *implicit* field" -- and closed ``GAP-CAD-002`` on that basis.
A B-rep kernel needs the other half, because a fillet, a hull/appendage seam, a deck edge
or a plate boundary is not the level set of anything: it is ``S1(u, v) = S2(s, t)``, two
NURBS patches that meet. Neither comes with a signed field, which is exactly why the
implicit route does not apply. ``GAP-CAD-A1-NURBS-INTERSECT`` records that gap and this
module is what closes it.

**The method.** The pair meets in a one-dimensional set, but the problem becomes an
ordinary root find the moment one patch is used as the parameter space:

    g(u, v) = (S1(u, v) - Q).n2,   Q = the closest point of S2 to S1(u, v),

i.e. the *oriented distance* from ``S1(u, v)`` to ``S2``. It vanishes on the intersection,
it is a signed length (so a residual like ``1e-9 L`` means what it says), and it is smooth
for transversal contact even while the closest point slides. That turns the curve into the
zero level set of a scalar field on a *flat rectangular parameter domain*, which marching
squares handles with no 3-D geometry: sample ``g`` on a grid, locate every sign change on a
grid edge, and chain those crossings cell to cell through edge keys -- so branch and loop
structure fall out of the cell topology rather than out of a proximity test (which is what
keeps two nearby branches apart).

The field is only the *seed* generator, not the answer. Every crossing is then polished on
the full system

    F(u, v, s, t) = S1(u, v) - S2(s, t) = 0

by a minimum-norm Gauss-Newton step on the 3 x 4 Jacobian ``[S1_u, S1_v, -S2_s, -S2_t]``.
That is what buys the accuracy: the march lands within ``O(h**2)`` of the curve, the polish
drives ``|F|`` to round-off, so the returned vertices lie on *both* surfaces to machine
precision rather than to grid precision. ``refine_levels`` then inserts parameter-space
midpoints and re-polishes, which halves the chord step per level and pulls the arc-length
error down as ``O(h**2)``; both the raw chord sum and the refined length are reported.

**Trim domains.** A closed intersection loop is the boundary of a trimmed face. For each
closed branch the module returns a :class:`TrimLoop` -- the same curve written in that
patch's own parameter space -- and a :class:`TrimDomain` that answers ``contains(u, v)`` for
the part of the patch lying on the inside of the *other* patch, using the sign of the same
oriented distance. Parameter area and surface area of that region come from masked
quadrature. The loop polygon also reports its exact shoelace area, but only under
``TrimLoop.area_reliable``: a loop that runs around a cylinder's seam, or that sits on a
constant-``v`` line of a lat/long sphere, encloses no parameter area and is flagged instead
of being handed a meaningless number.

**Honest failure modes.** Three are structural and reported, not hidden:

* *tangential contact* (a plane grazing a sphere, two spheres touching externally) produces
  no sign change, so the zero set is a point and not a curve. The march returns nothing and
  :attr:`NURBSIntersectionResult.tangential_candidate` reports that the two patches came
  within ``1e-6 L`` without crossing -- so a caller can tell "no curve" from "touched and I
  could not resolve it".
* *coincident patches* vanish on a whole region -- the intersection would be a surface. That
  raises :class:`NURBSIntersectionError` instead of returning a cloud of meaningless branches.
* *a feature narrower than one cell* is invisible. ``samples`` is the caller's dial and the
  grid actually used is echoed in ``NURBSIntersectionResult.diagnostics``.

Critical contact points, where the two normals are parallel (the crossing point of two
Viviani branches, osculating contact), make the 3 x 4 Jacobian rank-deficient; the polish
then falls back to a Tikhonov-damped step and *counts* those steps in
``diagnostics["damped_steps"]`` rather than swallowing them.
"""
from __future__ import annotations

import math
from collections import defaultdict
from dataclasses import dataclass, field, replace
from typing import Callable, Sequence

Vec2 = tuple[float, float]
Vec3 = tuple[float, float, float]

#: Grid *cells* per direction on the first patch. Raise until the answer stops moving.
DEFAULT_SAMPLES = (128, 40)
#: Midpoint insertions per segment after the march; each halves the chord step.
DEFAULT_REFINE_LEVELS = 2
#: Two branch ends closer than this times the patch scale are the same point.
DEFAULT_CLOSURE_TOLERANCE = 1e-9
#: A polished crossing further than this (times the scale) from ``F = 0`` is discarded.
DEFAULT_RESIDUAL_TOLERANCE = 1e-8
#: Tolerance of the closest-point Newton iteration, relative to the parameter span.
_PROJECTION_TOLERANCE = 1e-15
#: Gauss-Newton cap. Transversal contact converges in 3-5 steps; the rest is safety.
_POLISH_STEPS = 24
_CLOSEST_STEPS = 24
#: Coarse lattice used to seed the closest-point search on the second patch.
COARSE_SEEDS = (16, 8)
_TINY = 1e-300


class NURBSIntersectionError(ValueError):
    """Raised for input that cannot describe a NURBS x NURBS intersection."""


# --------------------------------------------------------------------------- #
# small vector helpers (the package stays standard-library only)
# --------------------------------------------------------------------------- #
def _sub(a: Vec3, b: Vec3) -> Vec3:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _scale(a: Vec3, s: float) -> Vec3:
    return (a[0] * s, a[1] * s, a[2] * s)


def _dot(a: Vec3, b: Vec3) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _cross(a: Vec3, b: Vec3) -> Vec3:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _norm(a: Vec3) -> float:
    return math.sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2])


def _chord(a: Vec3, b: Vec3) -> float:
    return _norm(_sub(a, b))


def _mid(a: Vec3, b: Vec3) -> Vec3:
    return (0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1]), 0.5 * (a[2] + b[2]))


# --------------------------------------------------------------------------- #
# surface adapter: what the march needs from a patch, and where it came from
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class SurfaceAdapter:
    """A patch normalised to ``domain`` + ``evaluate`` + ``jacobian`` + ``periodic``.

    ``derivative_source`` is ``"analytic"`` when the object exposed ``jacobian`` (every
    :class:`tensorcad.nurbs.NURBSSurface` does) and ``"finite-difference"`` when the
    derivatives had to be differenced. The distinction is kept because it bounds the
    accuracy of the Newton *steps*, even though it does not bound the accuracy of the
    *residual*: the residual is evaluated on the exact surface either way.
    """

    label: str
    raw_evaluate: Callable[[float, float], Vec3]
    raw_jacobian: Callable[[float, float], tuple[Vec3, Vec3]]
    u_range: Vec2
    v_range: Vec2
    periodic: tuple[bool, bool]
    derivative_source: str

    @property
    def u_span(self) -> float:
        return self.u_range[1] - self.u_range[0]

    @property
    def v_span(self) -> float:
        return self.v_range[1] - self.v_range[0]

    @staticmethod
    def _wrap_axis(value: float, span: Vec2, periodic: bool) -> float:
        if periodic:
            period = span[1] - span[0]
            return span[0] + math.fmod(math.fmod(value - span[0], period) + period, period)
        if value < span[0]:
            return span[0]
        if value > span[1]:
            return span[1]
        return value

    def wrap(self, u: float, v: float) -> Vec2:
        return (
            self._wrap_axis(u, self.u_range, self.periodic[0]),
            self._wrap_axis(v, self.v_range, self.periodic[1]),
        )

    def evaluate(self, u: float, v: float) -> Vec3:
        wu, wv = self.wrap(u, v)
        return self.raw_evaluate(wu, wv)

    def jacobian(self, u: float, v: float) -> tuple[Vec3, Vec3]:
        wu, wv = self.wrap(u, v)
        return self.raw_jacobian(wu, wv)

    def normal(self, u: float, v: float) -> Vec3 | None:
        """Unit normal, or ``None`` at a singular (zero-area) parameter point."""

        su, sv = self.jacobian(u, v)
        n = _cross(su, sv)
        length = _norm(n)
        if length <= _TINY:
            return None
        return _scale(n, 1.0 / length)

    def metric(self, u: float, v: float) -> float:
        """Area density ``|S_u x S_v|``: parameter area in, surface area out."""

        su, sv = self.jacobian(u, v)
        return _norm(_cross(su, sv))


def _adapt(surface: object, label: str) -> SurfaceAdapter:
    """Normalise a patch-like object; accepts anything with ``domain`` + ``evaluate``."""

    if not hasattr(surface, "evaluate"):
        raise NURBSIntersectionError(f"{label}: surface must expose evaluate(u, v)")
    ranges = getattr(surface, "domain", None)
    if ranges is None:
        raise NURBSIntersectionError(
            f"{label}: surface must expose a domain ((u0, u1), (v0, v1))"
        )
    if callable(ranges):
        ranges = ranges()
    try:
        (u0, u1), (v0, v1) = ranges
        u_range = (float(u0), float(u1))
        v_range = (float(v0), float(v1))
    except (TypeError, ValueError) as error:
        raise NURBSIntersectionError(
            f"{label}: surface.domain must be ((u0, u1), (v0, v1))"
        ) from error
    for low, high in (u_range, v_range):
        if not (math.isfinite(low) and math.isfinite(high)) or high <= low:
            raise NURBSIntersectionError(f"{label}: empty domain {(u_range, v_range)!r}")
    flags = tuple(getattr(surface, "periodic", (False, False)))
    if len(flags) != 2:
        raise NURBSIntersectionError(f"{label}: surface.periodic must be a pair of booleans")
    periodic = (bool(flags[0]), bool(flags[1]))

    raw = surface.evaluate  # type: ignore[attr-defined]

    def raw_evaluate(u: float, v: float) -> Vec3:
        point = raw(u, v)
        return (float(point[0]), float(point[1]), float(point[2]))

    supplied = getattr(surface, "jacobian", None)
    if callable(supplied):

        def raw_jacobian(u: float, v: float) -> tuple[Vec3, Vec3]:
            out = supplied(u, v)
            su, sv = out[0], out[1]
            return (
                (float(su[0]), float(su[1]), float(su[2])),
                (float(sv[0]), float(sv[1]), float(sv[2])),
            )

        source = "analytic"
    else:
        hu = 1e-6 * (u_range[1] - u_range[0])
        hv = 1e-6 * (v_range[1] - v_range[0])

        def raw_jacobian(u: float, v: float) -> tuple[Vec3, Vec3]:
            up, um = raw_evaluate(u + hu, v), raw_evaluate(u - hu, v)
            vp, vm = raw_evaluate(u, v + hv), raw_evaluate(u, v - hv)
            return (
                _scale(_sub(up, um), 0.5 / hu),
                _scale(_sub(vp, vm), 0.5 / hv),
            )

        source = "finite-difference"

    return SurfaceAdapter(
        label=label,
        raw_evaluate=raw_evaluate,
        raw_jacobian=raw_jacobian,
        u_range=u_range,
        v_range=v_range,
        periodic=periodic,
        derivative_source=source,
    )


# --------------------------------------------------------------------------- #
# closest point on a patch, and the oriented distance built on top of it
# --------------------------------------------------------------------------- #
def _closest_point(
    surface: SurfaceAdapter,
    point: Vec3,
    guess: Vec2,
    *,
    tolerance: float = _PROJECTION_TOLERANCE,
    steps: int = _CLOSEST_STEPS,
) -> tuple[Vec2, float]:
    """Gauss-Newton projection of ``point`` onto ``surface``, starting from ``guess``.

    Deliberately a *local* solve: the caller supplies the seed, because a global
    nearest-point query on a curved patch is a different and much harder problem, and the
    march only ever needs the branch it is already standing on.
    """

    u, v = surface.wrap(*guess)
    tu = tolerance * surface.u_span
    tv = tolerance * surface.v_span
    span_u, span_v = surface.u_span, surface.v_span
    for _ in range(steps):
        s = surface.evaluate(u, v)
        su, sv = surface.jacobian(u, v)
        residual = _sub(s, point)
        g1, g2 = _dot(su, residual), _dot(sv, residual)
        a11, a12, a22 = _dot(su, su), _dot(su, sv), _dot(sv, sv)
        det = a11 * a22 - a12 * a12
        trace = a11 + a22
        if trace <= _TINY:
            break
        if abs(det) <= 1e-12 * trace * trace:
            # A patch that is momentarily degenerate (a pole, a collapsed edge): damp
            # instead of dividing by a vanishing determinant.
            damp = 1e-8 * trace
            a11 += damp
            a22 += damp
            det = a11 * a22 - a12 * a12
            if abs(det) <= _TINY:
                break
        du = -(a22 * g1 - a12 * g2) / det
        dv = -(-a12 * g1 + a11 * g2) / det
        # A step larger than half the domain means the seed was on the wrong branch and
        # the iteration is jumping; clamp so it cannot wrap around a cylinder.
        du = max(-0.5 * span_u, min(0.5 * span_u, du))
        dv = max(-0.5 * span_v, min(0.5 * span_v, dv))
        u, v = surface.wrap(u + du, v + dv)
        if abs(du) <= tu and abs(dv) <= tv:
            break
    final = surface.evaluate(u, v)
    return (u, v), _chord(final, point)


def _seed_lattice(
    surface: SurfaceAdapter, counts: tuple[int, int]
) -> tuple[tuple[float, float, Vec3], ...]:
    nu, nv = max(2, int(counts[0])), max(2, int(counts[1]))
    u0, u1 = surface.u_range
    v0, v1 = surface.v_range
    seeds: list[tuple[float, float, Vec3]] = []
    for i in range(nu):
        u = u0 + (u1 - u0) * (i + 0.5) / nu
        for j in range(nv):
            v = v0 + (v1 - v0) * (j + 0.5) / nv
            seeds.append((u, v, surface.evaluate(u, v)))
    return tuple(seeds)


@dataclass
class _DistanceField:
    """Oriented distance from a point to ``target``: ``|g|`` is a length, the sign says
    which side of the target's tangent plane the point falls on.

    Seeding is *coherent*: the previous answer is tried first, then the nearest coarse
    lattice node. That keeps the march fast and keeps it on the branch it is standing on,
    which a global nearest-point search would not.
    """

    target: SurfaceAdapter
    seeds: tuple[tuple[float, float, Vec3], ...]
    steps: int = _CLOSEST_STEPS
    tolerance: float = _PROJECTION_TOLERANCE
    _last: Vec2 | None = None
    evaluations: int = 0

    def _nearest_seed(self, point: Vec3) -> Vec2:
        best_uv = (
            0.5 * (self.target.u_range[0] + self.target.u_range[1]),
            0.5 * (self.target.v_range[0] + self.target.v_range[1]),
        )
        best = math.inf
        for u, v, p in self.seeds:
            d = _chord(p, point)
            if d < best:
                best = d
                best_uv = (u, v)
        return best_uv

    def __call__(self, point: Vec3) -> tuple[float, Vec2]:
        candidates: list[Vec2] = []
        if self._last is not None:
            candidates.append(self._last)
        candidates.append(self._nearest_seed(point))
        best_uv = candidates[0]
        best_dist = math.inf
        for guess in candidates:
            uv, dist = _closest_point(
                self.target, point, guess, tolerance=self.tolerance, steps=self.steps
            )
            if dist < best_dist:
                best_dist, best_uv = dist, uv
            if best_dist <= 1e-13:
                break
        self._last = best_uv
        self.evaluations += 1
        foot = self.target.evaluate(*best_uv)
        normal = self.target.normal(*best_uv)
        if normal is None:
            # A singular point of the target (a pole): step one cell in, where the tangent
            # plane is actually defined.
            normal = self.target.normal(best_uv[0] + 1e-6 * self.target.u_span, best_uv[1])
        if normal is None:
            return best_dist, best_uv
        return _dot(_sub(point, foot), normal), best_uv

# --------------------------------------------------------------------------- #
# the 3 x 4 minimum-norm Gauss-Newton step
# --------------------------------------------------------------------------- #
def _det3(m: list[list[float]]) -> float:
    return (
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    )


def _solve3(m: list[list[float]], rhs: Vec3) -> Vec3 | None:
    """Solve the 3x3 system by Cramer's rule; ``None`` when it is singular."""

    det = _det3(m)
    if det == 0.0 or not math.isfinite(det):
        return None
    inv = 1.0 / det
    c00 = m[1][1] * m[2][2] - m[1][2] * m[2][1]
    c01 = -(m[0][1] * m[2][2] - m[0][2] * m[2][1])
    c02 = m[0][1] * m[1][2] - m[0][2] * m[1][1]
    c10 = -(m[1][0] * m[2][2] - m[1][2] * m[2][0])
    c11 = m[0][0] * m[2][2] - m[0][2] * m[2][0]
    c12 = -(m[0][0] * m[1][2] - m[0][2] * m[1][0])
    c20 = m[1][0] * m[2][1] - m[1][1] * m[2][0]
    c21 = -(m[0][0] * m[2][1] - m[0][1] * m[2][0])
    c22 = m[0][0] * m[1][1] - m[0][1] * m[1][0]
    return (
        inv * (c00 * rhs[0] + c01 * rhs[1] + c02 * rhs[2]),
        inv * (c10 * rhs[0] + c11 * rhs[1] + c12 * rhs[2]),
        inv * (c20 * rhs[0] + c21 * rhs[1] + c22 * rhs[2]),
    )


def _min_norm_step(
    columns: tuple[Vec3, Vec3, Vec3, Vec3], rhs: Vec3
) -> tuple[tuple[float, float, float, float], bool] | None:
    """``delta = -J^T (J J^T)^-1 F`` for the 3x4 Jacobian whose columns these are.

    Returns ``(delta, damped)``. ``damped`` is ``True`` when ``J J^T`` was numerically
    rank-deficient -- the two normals are parallel, i.e. the contact there is (nearly)
    tangential -- and a Tikhonov term had to be added to keep the step finite.
    """

    m = [[0.0] * 3 for _ in range(3)]
    for c in columns:
        for i in range(3):
            for j in range(3):
                m[i][j] += c[i] * c[j]
    trace = m[0][0] + m[1][1] + m[2][2]
    if trace <= _TINY:
        return None
    damped = False
    mean = trace / 3.0
    if abs(_det3(m)) <= 1e-10 * mean ** 3:
        damped = True
        lam = 1e-10 * mean
        for i in range(3):
            m[i][i] += lam
        if abs(_det3(m)) <= _TINY:
            return None
    x = _solve3(m, rhs)
    if x is None:
        return None
    delta = (
        -_dot(columns[0], x),
        -_dot(columns[1], x),
        -_dot(columns[2], x),
        -_dot(columns[3], x),
    )
    return delta, damped


def _polish(
    surface_a: SurfaceAdapter,
    surface_b: SurfaceAdapter,
    u: float,
    v: float,
    s: float,
    t: float,
    *,
    scale: float,
    steps: int = _POLISH_STEPS,
) -> tuple[Vec2, Vec2, float, int]:
    """Drive ``S1(u, v) - S2(s, t)`` to zero, keeping the step inside the tangent band.

    Returns ``(params_a, params_b, best_gap, damped_steps)``. ``best_gap`` is the *smallest*
    ``|S1 - S2|`` seen, not the last one: a damped step can chatter at round-off and the
    best iterate is the honest answer.
    """

    best = (u, v, s, t, math.inf)
    damped_steps = 0
    limits = (
        0.25 * surface_a.u_span,
        0.25 * surface_a.v_span,
        0.25 * surface_b.u_span,
        0.25 * surface_b.v_span,
    )
    target = 1e-15 * scale
    for _ in range(steps):
        gap = _chord(surface_a.evaluate(u, v), surface_b.evaluate(s, t))
        if gap < best[4]:
            best = (u, v, s, t, gap)
        if gap <= target:
            break
        au, av = surface_a.jacobian(u, v)
        bs, bt = surface_b.jacobian(s, t)
        step = _min_norm_step((au, av, _scale(bs, -1.0), _scale(bt, -1.0)),
                              _sub(surface_a.evaluate(u, v), surface_b.evaluate(s, t)))
        if step is None:
            break
        delta, damped = step
        if damped:
            damped_steps += 1
        du, dv, ds, dt = (min(limit, max(-limit, d)) for d, limit in zip(delta, limits))
        u, v = surface_a.wrap(u + du, v + dv)
        s, t = surface_b.wrap(s + ds, t + dt)
    return (best[0], best[1]), (best[2], best[3]), best[4], damped_steps


# --------------------------------------------------------------------------- #
# the march: grid -> crossings -> chained branches
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class _Crossing:
    """One sign change on a grid edge, already polished onto both surfaces."""

    key: tuple[str, int, int]
    params_a: Vec2
    params_b: Vec2
    point: Vec3
    gap: float
    damped: int


def _lerp_params(
    surface: SurfaceAdapter, p0: Vec2, p1: Vec2, weight: float
) -> Vec2:
    """Linear interpolation of a parameter pair, unwrapping periodic axes the short way."""

    out: list[float] = []
    for axis, (low, high) in enumerate((surface.u_range, surface.v_range)):
        a, b = p0[axis], p1[axis]
        if surface.periodic[axis]:
            period = high - low
            delta = b - a
            delta -= period * round(delta / period)
            b = a + delta
        out.append(a + weight * (b - a))
    return (out[0], out[1])


def _mid_params(surface: SurfaceAdapter, p0: Vec2, p1: Vec2) -> Vec2:
    return _lerp_params(surface, p0, p1, 0.5)


def _characteristic_scale(*surfaces: SurfaceAdapter) -> float:
    """Diameter of the union of the patches' four-corner boxes -- every tolerance unit."""

    corners: list[Vec3] = []
    for surface in surfaces:
        for u in (surface.u_range[0], surface.u_range[1]):
            for v in (surface.v_range[0], surface.v_range[1]):
                corners.append(surface.evaluate(u, v))
    centre = (
        sum(c[0] for c in corners) / len(corners),
        sum(c[1] for c in corners) / len(corners),
        sum(c[2] for c in corners) / len(corners),
    )
    scale = max((_chord(c, centre) for c in corners), default=0.0) * 2.0
    if not math.isfinite(scale) or scale <= 0.0:
        raise NURBSIntersectionError("the two patches have zero extent; nothing to cut")
    return scale


def _grid_values(
    surface: SurfaceAdapter,
    field: _DistanceField,
    cells: tuple[int, int],
) -> tuple[list[list[float]], list[float], list[float]]:
    """Sample the oriented distance on the node grid, plus the node parameter values."""

    u_cells, v_cells = cells
    u_nodes = u_cells if surface.periodic[0] else u_cells + 1
    v_nodes = v_cells if surface.periodic[1] else v_cells + 1
    du = surface.u_span / u_cells
    dv = surface.v_span / v_cells
    u_values = [surface.u_range[0] + i * du for i in range(u_nodes)]
    v_values = [surface.v_range[0] + j * dv for j in range(v_nodes)]
    values = [
        [field(surface.evaluate(u_values[i], v_values[j]))[0] for j in range(v_nodes)]
        for i in range(u_nodes)
    ]
    return values, u_values, v_values


def _march(
    surface_a: SurfaceAdapter,
    surface_b: SurfaceAdapter,
    field: _DistanceField,
    cells: tuple[int, int],
    values: list[list[float]],
    u_values: list[float],
    v_values: list[float],
    *,
    scale: float,
    residual_tolerance: float,
) -> tuple[list[tuple[list[_Crossing], bool]], dict[str, float]]:
    """Marching squares on the parameter grid, then polish every crossing."""

    u_cells, v_cells = cells
    u_nodes, v_nodes = len(u_values), len(v_values)
    du = surface_a.u_span / u_cells
    dv = surface_a.v_span / v_cells
    periodic_u, periodic_v = surface_a.periodic

    def iu(i: int) -> int:
        return i % u_nodes if periodic_u else i

    def jv(j: int) -> int:
        return j % v_nodes if periodic_v else j

    def sign(value: float) -> int:
        return 1 if value >= 0.0 else -1

    crossings: dict[tuple[str, int, int], _Crossing] = {}
    dropped = 0
    damped_total = 0

    def crossing(key: tuple[str, int, int], p0: Vec2, g0: float, p1: Vec2, g1: float) -> _Crossing | None:
        nonlocal dropped, damped_total
        if key in crossings:
            return crossings[key]
        weight = g0 / (g0 - g1) if g0 != g1 else 0.5
        weight = min(1.0, max(0.0, weight))
        guess_a = _lerp_params(surface_a, p0, p1, weight)
        point = surface_a.evaluate(*guess_a)
        _, seed_b = field(point)
        params_a, params_b, gap, damped = _polish(
            surface_a, surface_b, guess_a[0], guess_a[1], seed_b[0], seed_b[1],
            scale=scale,
        )
        damped_total += damped
        if not math.isfinite(gap) or gap > residual_tolerance * scale:
            dropped += 1
            return None
        record = _Crossing(
            key=key,
            params_a=params_a,
            params_b=params_b,
            point=surface_a.evaluate(*params_a),
            gap=gap,
            damped=damped,
        )
        crossings[key] = record
        return record

    adjacency: dict[tuple[str, int, int], list[tuple[str, int, int]]] = defaultdict(list)

    for i in range(u_cells):
        i1 = iu(i + 1)
        for j in range(v_cells):
            j1 = jv(j + 1)
            corners = ((i, j), (i1, j), (i1, j1), (i, j1))
            edge_keys = (("u", i, j), ("v", i1, j), ("u", i, j1), ("v", i, j))
            present: list[tuple[str, int, int]] = []
            for k in range(4):
                (ia, ja), (ib, jb) = corners[k], corners[(k + 1) % 4]
                ga, gb = values[ia][ja], values[ib][jb]
                if sign(ga) == sign(gb):
                    continue
                pa = (u_values[ia], v_values[ja])
                pb = (u_values[ib], v_values[jb])
                record = crossing(edge_keys[k], pa, ga, pb, gb)
                if record is not None:
                    present.append(edge_keys[k])
            if len(present) == 2:
                adjacency[present[0]].append(present[1])
                adjacency[present[1]].append(present[0])
            elif len(present) == 4:
                # A saddle: the cell centre decides which pairs belong together, i.e. the
                # isolated corners are the ones whose sign differs from the centre.
                centre_value = field(
                    surface_a.evaluate(u_values[i] + 0.5 * du, v_values[j] + 0.5 * dv)
                )[0]
                centre_sign = sign(centre_value)
                pairs: list[tuple[tuple[str, int, int], tuple[str, int, int]]] = []
                for k in range(4):
                    ia, ja = corners[k]
                    if sign(values[ia][ja]) == centre_sign:
                        continue
                    left, right = edge_keys[(k - 1) % 4], edge_keys[k]
                    if left in present and right in present:
                        pairs.append((left, right))
                for left, right in pairs:
                    adjacency[left].append(right)
                    adjacency[right].append(left)

    chains: list[tuple[list[_Crossing], bool]] = []
    visited: set[tuple[str, int, int]] = set()

    def walk(start: tuple[str, int, int]) -> tuple[list[tuple[str, int, int]], bool]:
        chain = [start]
        visited.add(start)
        previous, current = None, start
        while True:
            options = [k for k in adjacency.get(current, ()) if k != previous]
            if not options:
                return chain, False
            nxt = options[0]
            if nxt == start:
                return chain, len(chain) > 2
            if nxt in visited:
                return chain, False
            chain.append(nxt)
            visited.add(nxt)
            previous, current = current, nxt

    for key in list(adjacency):
        if key in visited or len(adjacency[key]) != 1:
            continue
        keys, _ = walk(key)
        chains.append(([crossings[k] for k in keys], False))
    for key in list(adjacency):
        if key in visited or len(adjacency[key]) < 2:
            continue
        keys, loop_closed = walk(key)
        chains.append(([crossings[k] for k in keys], loop_closed))

    diagnostics = {
        "grid_cells": float(u_cells),
        "grid_cells_v": float(v_cells),
        "crossings": float(len(crossings)),
        "dropped_crossings": float(dropped),
        "damped_steps": float(damped_total),
        "min_abs_field": min(abs(g) for row in values for g in row),
    }
    return chains, diagnostics


# --------------------------------------------------------------------------- #
# refinement: halve the chord step, re-polish, halve again
# --------------------------------------------------------------------------- #
def _refine(
    surface_a: SurfaceAdapter,
    surface_b: SurfaceAdapter,
    records: Sequence[_Crossing],
    *,
    closed_by: str,
    scale: float,
    levels: int,
) -> tuple[list[_Crossing], int]:
    """Insert parameter-space midpoints and project each onto ``S1 = S2``.

    The chord error of a segment of length ``h`` on a curve of curvature ``k`` is
    ``k**2 h**3 / 24``, so the *relative* length error behaves like ``h**2`` and each
    level should cut it by four. That is a prediction the example measures, not a claim
    this docstring gets to make on its own.

    The closing segment of a *positional* (seam) closure is deliberately **not** refined:
    its two parameters sit at opposite ends of the domain, so their midpoint is a totally
    different point of the surface. Only its chord is added to the length.
    """

    current = list(records)
    damped_total = 0
    for _ in range(levels):
        refined: list[_Crossing] = []
        for k in range(len(current) - 1):
            left, right = current[k], current[k + 1]
            refined.append(left)
            mid_a = _mid_params(surface_a, left.params_a, right.params_a)
            mid_b = _mid_params(surface_b, left.params_b, right.params_b)
            params_a, params_b, gap, damped = _polish(
                surface_a, surface_b, mid_a[0], mid_a[1], mid_b[0], mid_b[1], scale=scale
            )
            damped_total += damped
            if not math.isfinite(gap) or gap > 1e-6 * scale:
                # The midpoint could not be projected; leaving it out keeps the segment
                # honest (one long chord) instead of bending the curve onto a wrong branch.
                continue
            refined.append(
                _Crossing(
                    key=left.key,
                    params_a=params_a,
                    params_b=params_b,
                    point=surface_a.evaluate(*params_a),
                    gap=gap,
                    damped=damped,
                )
            )
        refined.append(current[-1])
        if closed_by == "loop" and len(current) > 2:
            left, right = current[-1], current[0]
            mid_a = _mid_params(surface_a, left.params_a, right.params_a)
            mid_b = _mid_params(surface_b, left.params_b, right.params_b)
            params_a, params_b, gap, damped = _polish(
                surface_a, surface_b, mid_a[0], mid_a[1], mid_b[0], mid_b[1], scale=scale
            )
            damped_total += damped
            if math.isfinite(gap) and gap <= 1e-6 * scale:
                refined.append(
                    _Crossing(
                        key=left.key,
                        params_a=params_a,
                        params_b=params_b,
                        point=surface_a.evaluate(*params_a),
                        gap=gap,
                        damped=damped,
                    )
                )
        current = refined
    return current, damped_total


# --------------------------------------------------------------------------- #
# results
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class NURBSIntersectionCurve:
    """One branch of ``S1 = S2``: the points, both parameter traces, and how it closed."""

    index: int
    points: tuple[Vec3, ...]
    params_a: tuple[Vec2, ...]
    params_b: tuple[Vec2, ...]
    closed: bool
    closed_by: str
    length: float
    chord_length: float
    max_gap: float

    def __post_init__(self) -> None:
        if len(self.points) != len(self.params_a) or len(self.points) != len(self.params_b):
            raise NURBSIntersectionError("points and both parameter traces must agree in length")

    @property
    def closure_defect(self) -> float:
        """Distance between the two ends: zero for a branch the march closed up."""

        if len(self.points) < 2:
            return 0.0
        return _chord(self.points[0], self.points[-1])


@dataclass(frozen=True)
class TrimLoop:
    """The intersection curve written in one patch's own parameter space.

    ``area_uv`` is the shoelace area of the traced polygon and is only meaningful when
    ``area_reliable`` -- a loop that closes across a periodic seam, or that lies along a
    constant-``v`` line, spans the domain instead of enclosing part of it.
    """

    params: tuple[Vec2, ...]
    closed: bool
    area_uv: float
    area_reliable: bool
    seam_spanning: bool
    perimeter_uv: float


@dataclass(frozen=True)
class TrimDomain:
    """The part of one patch that lies on the inside of the other -- i.e. what a trim cuts.

    ``contains`` is answered by the sign of the oriented distance to the other patch (the
    same field the march uses), not by a polygon test, so it stays correct for loops that
    run around a seam. The loop polygons are kept for callers that need the boundary
    itself (B-rep edge construction, drawing, area bookkeeping).
    """

    label: str
    loops: tuple[TrimLoop, ...]
    domain: tuple[Vec2, Vec2]
    periodic: tuple[bool, bool]
    convention: str
    _inside: Callable[[float, float], bool] = field(repr=False, compare=False)
    _metric: Callable[[float, float], float] = field(repr=False, compare=False)

    def contains(self, u: float, v: float) -> bool:
        """Is the patch point ``(u, v)`` inside the trimmed (kept) region?"""

        return bool(self._inside(float(u), float(v)))

    def _masked_integral(
        self, samples: Sequence[int], *, density: bool, refine: int
    ) -> float:
        """Sum a cell-wise integrand over the kept region, bisecting the boundary cells.

        A plain midpoint rule loses ``O(h)`` on the boundary of the region, because the
        boundary is a curve and not an area. Splitting only the cells whose four corners
        disagree (there are ``O(h)`` of them) and sampling those ``4**refine`` ways removes
        most of that bias at a small cost, and -- more to the point -- it makes
        ``area_uv`` a number whose accuracy can be *checked* against an analytic value
        instead of assumed. It does not remove the bias entirely: the residual is
        first-order on the boundary, which is why the reported areas are good to a few
        1e-3 relative, not to round-off.
        """

        u_cells, v_cells = int(samples[0]), int(samples[1])
        u0, u1 = self.domain[0]
        v0, v1 = self.domain[1]
        du = (u1 - u0) / u_cells
        dv = (v1 - v0) / v_cells
        sub = 2 ** max(0, int(refine))
        total = 0.0

        def density_at(u: float, v: float) -> float:
            return self._metric(u, v) if density else 1.0

        for i in range(u_cells):
            ua = u0 + i * du
            for j in range(v_cells):
                va = v0 + j * dv
                corners = (
                    self._inside(ua, va),
                    self._inside(ua + du, va),
                    self._inside(ua + du, va + dv),
                    self._inside(ua, va + dv),
                )
                if all(corners) or not any(corners):
                    if corners[0]:
                        total += density_at(ua + 0.5 * du, va + 0.5 * dv) * du * dv
                    continue
                su = du / sub
                sv = dv / sub
                for a in range(sub):
                    for b in range(sub):
                        uu = ua + (a + 0.5) * su
                        vv = va + (b + 0.5) * sv
                        if self._inside(uu, vv):
                            total += density_at(uu, vv) * su * sv
        return total

    def area_uv(
        self, samples: Sequence[int] = (64, 32), *, refine: int = 2
    ) -> float:
        """Parameter area of the trimmed region (masked, boundary-refined quadrature)."""

        return self._masked_integral(samples, density=False, refine=refine)

    def surface_area(
        self, samples: Sequence[int] = (64, 32), *, refine: int = 2
    ) -> float:
        """Surface area of the trimmed region, ``int int |S_u x S_v| du dv`` over it."""

        return self._masked_integral(samples, density=True, refine=refine)


@dataclass(frozen=True)
class NURBSIntersectionResult:
    """Everything the march found on the pair, plus what it cost and how it failed."""

    curves: tuple[NURBSIntersectionCurve, ...]
    trim_a: TrimDomain
    trim_b: TrimDomain
    samples: tuple[int, int]
    refine_levels: int
    scale: float
    diagnostics: dict[str, float]

    @property
    def total_length(self) -> float:
        return sum(curve.length for curve in self.curves)

    @property
    def total_chord_length(self) -> float:
        return sum(curve.chord_length for curve in self.curves)

    @property
    def closed_curves(self) -> tuple[NURBSIntersectionCurve, ...]:
        return tuple(curve for curve in self.curves if curve.closed)

    @property
    def open_curves(self) -> tuple[NURBSIntersectionCurve, ...]:
        return tuple(curve for curve in self.curves if not curve.closed)

    @property
    def point_count(self) -> int:
        return sum(len(curve.points) for curve in self.curves)

    @property
    def max_gap(self) -> float:
        return max((curve.max_gap for curve in self.curves), default=0.0)

    @property
    def tangential_candidate(self) -> bool:
        """True when the patches came within ``1e-6 L`` but no curve came back.

        That is a touch, not a curve: the march has nothing to chain, and the honest
        answer is "one point, not representable as a branch". A caller that needs to know
        whether ``curves == ()`` means "disjoint" or "touching" reads this flag.
        """

        return bool(
            not self.curves
            and self.diagnostics.get("min_abs_field", math.inf) <= 1e-6 * self.scale
        )


def _trace_loop(
    params: Sequence[Vec2], closed: bool, surface: SurfaceAdapter
) -> TrimLoop:
    span_u, span_v = surface.u_span, surface.v_span
    pairs = list(zip(params, list(params)[1:] + ([params[0]] if closed and params else [])))
    seam = any(abs(b[0] - a[0]) > 0.5 * span_u or abs(b[1] - a[1]) > 0.5 * span_v
               for a, b in pairs)
    area = 0.0
    if len(params) >= 3:
        for (u0, v0), (u1, v1) in pairs:
            area += u0 * v1 - u1 * v0
        area = abs(0.5 * area)
    perimeter = sum(math.hypot(b[0] - a[0], b[1] - a[1]) for a, b in pairs)
    floor = 1e-9 * span_u * span_v
    return TrimLoop(
        params=tuple(params),
        closed=closed,
        area_uv=area,
        area_reliable=bool(closed and not seam and area > floor),
        seam_spanning=bool(seam),
        perimeter_uv=perimeter,
    )


def _trim_domain(
    surface: SurfaceAdapter,
    curves: Sequence[NURBSIntersectionCurve],
    *,
    use_first: bool,
    field_fn: _DistanceField,
) -> TrimDomain:
    """Loop traces plus a membership test, for the side of ``surface`` inside the other.

    ``field_fn`` must be the oriented distance *to the other patch*, so ``< 0`` means
    "inside". Passing the field of the wrong patch would return a trim domain that is the
    complement of the right one, which is why the two calls below pass different fields.
    """

    if use_first:
        loops = [_trace_loop(curve.params_a, curve.closed, surface) for curve in curves if curve.closed]
    else:
        loops = [_trace_loop(curve.params_b, curve.closed, surface) for curve in curves if curve.closed]

    def inside(u: float, v: float) -> bool:
        return field_fn(surface.evaluate(u, v))[0] < 0.0

    return TrimDomain(
        label=surface.label,
        loops=tuple(loops),
        domain=(surface.u_range, surface.v_range),
        periodic=surface.periodic,
        convention="inside = oriented distance to the other patch is negative",
        _inside=inside,
        _metric=surface.metric,
    )


# --------------------------------------------------------------------------- #
# the public entry point
# --------------------------------------------------------------------------- #
def nurbs_surface_intersection(
    surface_a: object,
    surface_b: object,
    *,
    samples: Sequence[int] = DEFAULT_SAMPLES,
    refine_levels: int = DEFAULT_REFINE_LEVELS,
    closure_tolerance: float = DEFAULT_CLOSURE_TOLERANCE,
    residual_tolerance: float = DEFAULT_RESIDUAL_TOLERANCE,
) -> NURBSIntersectionResult:
    """Intersect two parametric patches and trim each one by the other.

    ``surface_a`` and ``surface_b`` each need ``domain`` (``((u0, u1), (v0, v1))``) and
    ``evaluate(u, v)``; a ``jacobian(u, v) -> (S_u, S_v)`` is used when present (every
    :class:`tensorcad.nurbs.NURBSSurface` has it) and differenced when not. ``periodic``
    (a pair of booleans) is honoured when declared -- an angle parameter whose seam would
    otherwise tear a closed loop into an open chain.

    ``samples`` is the number of grid *cells* in ``u`` and ``v`` on ``surface_a``; the
    oriented-distance field to ``surface_b`` is sampled there. ``refine_levels`` is how
    many times midpoints are inserted and re-projected (0-4). ``closure_tolerance`` and
    ``residual_tolerance`` are relative to the patches' own scale ``L`` (their bounding
    diameter), so nothing here is in absolute units.

    Raises :class:`NURBSIntersectionError` for a degenerate request (no domain, zero
    extent) and when the two patches coincide -- the latter because the intersection is
    then a *surface*, and returning a cloud of branches would be a lie.
    """

    if len(samples) != 2:
        raise NURBSIntersectionError("samples must be a pair of cell counts")
    u_cells, v_cells = int(samples[0]), int(samples[1])
    if u_cells < 8 or v_cells < 4:
        raise NURBSIntersectionError(
            f"samples must be at least (8, 4), got {(u_cells, v_cells)!r}"
        )
    if not 0 <= refine_levels <= 4:
        raise NURBSIntersectionError(f"refine_levels must be between 0 and 4, got {refine_levels!r}")
    for name, value in (("closure_tolerance", closure_tolerance), ("residual_tolerance", residual_tolerance)):
        if not math.isfinite(value) or value <= 0.0:
            raise NURBSIntersectionError(f"{name} must be positive, got {value!r}")

    adapter_a = _adapt(surface_a, "surface_a")
    adapter_b = _adapt(surface_b, "surface_b")
    scale = _characteristic_scale(adapter_a, adapter_b)

    field_b = _DistanceField(adapter_b, _seed_lattice(adapter_b, COARSE_SEEDS))
    values, u_values, v_values = _grid_values(adapter_a, field_b, (u_cells, v_cells))

    flat = [g for row in values for g in row]
    min_abs_field = min(abs(g) for g in flat)
    coincident_fraction = sum(1 for g in flat if abs(g) <= 1e-9 * scale) / len(flat)
    if coincident_fraction > 0.25:
        raise NURBSIntersectionError(
            "the two patches coincide over "
            f"{coincident_fraction:.0%} of the sampled domain: the intersection is a "
            "surface, not a curve"
        )

    sign_changes = 0
    for i in range(len(u_values) - (0 if adapter_a.periodic[0] else 1)):
        i1 = (i + 1) % len(u_values)
        for j in range(len(v_values) - (0 if adapter_a.periodic[1] else 1)):
            j1 = (j + 1) % len(v_values)
            for ga, gb in (
                (values[i][j], values[i1][j]),
                (values[i1][j], values[i1][j1]),
                (values[i1][j1], values[i][j1]),
                (values[i][j1], values[i][j]),
            ):
                if (ga >= 0.0) != (gb >= 0.0):
                    sign_changes += 1

    chains, march_diagnostics = _march(
        adapter_a, adapter_b, field_b, (u_cells, v_cells), values, u_values, v_values,
        scale=scale, residual_tolerance=residual_tolerance,
    )

    curves: list[NURBSIntersectionCurve] = []
    damped_total = march_diagnostics["damped_steps"]
    for records, loop_closed in chains:
        if len(records) < 2:
            continue
        points = [record.point for record in records]
        closed_by = ""
        if loop_closed:
            closed_by = "loop"
        elif _chord(points[0], points[-1]) <= closure_tolerance * scale:
            closed_by = "seam"
        closed = closed_by != ""
        refined, damped = _refine(
            adapter_a, adapter_b, records,
            closed_by=closed_by, scale=scale, levels=refine_levels,
        )
        damped_total += damped
        refined_points = [record.point for record in refined]
        chord_length = sum(
            _chord(a, b) for a, b in zip(refined_points, refined_points[1:])
        )
        length = chord_length
        if closed:
            length += _chord(refined_points[-1], refined_points[0])
        curves.append(
            NURBSIntersectionCurve(
                index=len(curves),
                points=tuple(refined_points),
                params_a=tuple(record.params_a for record in refined),
                params_b=tuple(record.params_b for record in refined),
                closed=closed,
                closed_by=closed_by,
                length=length,
                chord_length=chord_length,
                max_gap=max(record.gap for record in refined),
            )
        )

    # Membership is a sign question, not a precision one, so the trim fields run a cheap
    # projection: same surface, looser iteration, and the raster is thousands of calls.
    fast_b = _DistanceField(
        adapter_b, field_b.seeds, steps=4, tolerance=1e-8
    )
    fast_a = _DistanceField(
        adapter_a, _seed_lattice(adapter_a, COARSE_SEEDS), steps=4, tolerance=1e-8
    )
    trim_a = _trim_domain(adapter_a, curves, use_first=True, field_fn=fast_b)
    trim_b = _trim_domain(adapter_b, curves, use_first=False, field_fn=fast_a)

    # Tangential contact leaves a cloud of crossings that all polish onto the same point:
    # a trace of round-off, not a curve. Dropping it (and counting it) is what lets
    # ``tangential_candidate`` answer "touched" instead of returning a zero-length branch.
    collapsed = [curve for curve in curves if curve.length <= 1e-9 * scale]
    if collapsed:
        curves = [curve for curve in curves if curve.length > 1e-9 * scale]
        for position, curve in enumerate(curves):
            curves[position] = replace(curve, index=position)

    diagnostics = {
        "scale": scale,
        "collapsed_branches": float(len(collapsed)),
        "grid_cells": float(u_cells),
        "grid_cells_v": float(v_cells),
        "field_evaluations": float(field_b.evaluations),
        "sign_changes": float(sign_changes),
        "crossings": march_diagnostics["crossings"],
        "dropped_crossings": march_diagnostics["dropped_crossings"],
        "damped_steps": float(damped_total),
        "min_abs_field": min_abs_field,
        "coincident_fraction": coincident_fraction,
        "closure_tolerance": closure_tolerance,
        "residual_tolerance": residual_tolerance,
        "branches": float(len(curves)),
    }
    return NURBSIntersectionResult(
        curves=tuple(curves),
        trim_a=trim_a,
        trim_b=trim_b,
        samples=(u_cells, v_cells),
        refine_levels=refine_levels,
        scale=scale,
        diagnostics=diagnostics,
    )


__all__ = [
    "DEFAULT_CLOSURE_TOLERANCE",
    "DEFAULT_REFINE_LEVELS",
    "DEFAULT_RESIDUAL_TOLERANCE",
    "DEFAULT_SAMPLES",
    "NURBSIntersectionCurve",
    "NURBSIntersectionError",
    "NURBSIntersectionResult",
    "SurfaceAdapter",
    "TrimDomain",
    "TrimLoop",
    "nurbs_surface_intersection",
]
