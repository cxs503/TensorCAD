"""Dependency-free IGES reader: B-spline surfaces/curves, trimming, tessellation.

Why this module exists
---------------------
The baseline hull fleet ships as IGES files published by the workshops that
built the hulls (Tokyo 2015 / SIMMAN 2020).  Those files are the *only*
authoritative geometry for KCS, KVLCC2, JBC and ONRT, but the repository is
standard-library only by charter (``pyproject.toml``: ``dependencies = []``),
so an OpenCASCADE bridge would break the package contract.  This module reads
the small, well-defined slice of IGES 5.3 that those files actually use and
turns it into point samples; the meshing itself lives in
:mod:`tensorcad.fleet_baseline`.

Supported entities
------------------
=================  ============================================================
102                composite curve (used only for trim bookkeeping)
110                line
126                rational B-spline curve (3D or (u,v) parameter space)
128                rational B-spline surface
141                boundary (list of model-space / parameter-space curve pairs)
142                curve on a parametric surface
143                bounded surface (basis surface + one boundary)
144                trimmed parametric surface (basis surface + outer boundary)
=================  ============================================================

Deliberate limits
-----------------
* Only the delimiters ``','`` ``'.'`` ``';'`` are handled (every workshop file
  checked in this repository uses them); anything else raises :class:`IgesError`
  rather than silently mis-parsing.
* No transformation-matrix, view or colour processing: the hull files are
  already in model space, metres, origin at the midship section.
* Trim loops that are not representable as a closed chain of curve-on-surface
  segments are reported, not repaired.

Nothing here is a CAD kernel: it evaluates published geometry, it does not
build or heal surfaces.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Iterable, Sequence

__all__ = [
    "IgesError",
    "DirectoryEntry",
    "BSplineCurve",
    "BSplineSurface",
    "Boundary",
    "CurveOnSurface",
    "BoundedSurface",
    "IgesFile",
    "read_iges",
    "ENTITY_LINE",
    "ENTITY_COMPOSITE_CURVE",
    "ENTITY_BSPLINE_CURVE",
    "ENTITY_BSPLINE_SURFACE",
    "ENTITY_BOUNDARY",
    "ENTITY_CURVE_ON_SURFACE",
    "ENTITY_BOUNDED_SURFACE",
    "ENTITY_TRIMMED_SURFACE",
]

ENTITY_LINE = 110
ENTITY_COMPOSITE_CURVE = 102
ENTITY_BSPLINE_CURVE = 126
ENTITY_BSPLINE_SURFACE = 128
ENTITY_BOUNDARY = 141
ENTITY_CURVE_ON_SURFACE = 142
ENTITY_BOUNDED_SURFACE = 143
ENTITY_TRIMMED_SURFACE = 144

#: Type 126 is the **rational B-spline curve**, not a surface: IGES 5.3 defines no
#: non-rational surface entity, and every 126 record in the workshop archives
#: (KCS 0, JBC 172, KVLCC2 304, ONRT 820) decodes here as a 1-D-parameter curve
#: with a 3-D control net -- ``BSplineCurve.from_parameters`` has done so since the
#: first reader revision.  An earlier revision of this module declared a constant
#: ``ENTITY_BSPLINE_SURFACE_POLYNOMIAL = 126`` and reported these curves as
#: "unreadable surfaces"; the name was wrong twice (126 is a curve, and it is
#: decoded) and the numbers it produced were attributed to a reader gap.  See
#: ``docs/baseline_fleet_type126_falsified.md`` for the measurement that settled it.
ENTITY_BSPLINE_CURVE_ALIAS = ENTITY_BSPLINE_CURVE

#: Section letters of an IGES 5.x file; column 73 (0-based index 72).
_SECTIONS = "SGDPT"


class IgesError(ValueError):
    """Raised when a file cannot be read as IGES 5.3 without guessing."""


# --------------------------------------------------------------------------- #
# tokenising
# --------------------------------------------------------------------------- #
def _strip_number(token: str) -> float:
    """IGES writes exponent form as ``1.0D+3``; ``float`` needs an ``E``."""
    text = token.strip()
    if not text:
        raise IgesError("empty numeric field")
    if "D" in text or "d" in text:
        text = text.replace("D", "E").replace("d", "e")
    try:
        return float(text)
    except ValueError as exc:  # pragma: no cover - malformed input
        raise IgesError(f"not a number: {token!r}") from exc


def _tokens(payload: str) -> list[str]:
    return [t.strip() for t in payload.split(",") if t.strip()]


# --------------------------------------------------------------------------- #
# B-spline evaluation (de Boor, rational form)
# --------------------------------------------------------------------------- #
def _find_span(n: int, degree: int, knots: Sequence[float], value: float) -> int:
    """Index of the knot span containing ``value`` (The NURBS Book, A2.1)."""
    lo, hi = degree, n  # n = count-1
    if value >= knots[hi]:
        return hi - 1
    if value <= knots[lo]:
        return lo
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if value < knots[mid]:
            hi = mid
        else:
            lo = mid
    return lo


def _basis(n: int, degree: int, knots: Sequence[float], value: float) -> tuple[list[float], int]:
    """Non-zero basis functions ``N[span-degree .. span]`` at ``value``."""
    span = _find_span(n, degree, knots, value)
    values = [0.0] * (degree + 1)
    left = [0.0] * (degree + 1)
    right = [0.0] * (degree + 1)
    values[0] = 1.0
    for j in range(1, degree + 1):
        left[j] = value - knots[span + 1 - j]
        right[j] = knots[span + j] - value
        saved = 0.0
        for r in range(j):
            denominator = right[r + 1] + left[j - r]
            temp = 0.0 if denominator == 0.0 else values[r] / denominator
            values[r] = saved + right[r + 1] * temp
            saved = left[j - r] * temp
        values[j] = saved
    return values, span


@dataclass(frozen=True)
class BSplineCurve:
    """Rational B-spline curve of arbitrary dimension (2 for (u,v), 3 for model)."""

    degree: int
    knots: tuple[float, ...]
    control: tuple[tuple[float, ...], ...]
    weights: tuple[float, ...]
    dim: int

    @property
    def domain(self) -> tuple[float, float]:
        return self.knots[self.degree], self.knots[len(self.control)]

    def evaluate(self, t: float) -> tuple[float, ...]:
        count = len(self.control)
        value = min(max(t, self.knots[self.degree]), self.knots[count])
        basis, span = _basis(count - 1, self.degree, self.knots, value)
        numerator = [0.0] * self.dim
        denominator = 0.0
        for index in range(self.degree + 1):
            weight = self.weights[span - self.degree + index] * basis[index]
            point = self.control[span - self.degree + index]
            denominator += weight
            for axis in range(self.dim):
                numerator[axis] += weight * point[axis]
        if denominator == 0.0:
            return tuple(numerator[axis] for axis in range(self.dim))
        return tuple(numerator[axis] / denominator for axis in range(self.dim))

    def sample(self, count: int) -> list[tuple[float, ...]]:
        lo, hi = self.domain
        if count < 2:
            raise ValueError("need at least two samples")
        return [self.evaluate(lo + (hi - lo) * i / (count - 1)) for i in range(count)]

    @classmethod
    def from_parameters(cls, params: Sequence[float]) -> "BSplineCurve":
        """Build from the raw type-126 parameter list (entity type included).

        Field order is IGES 5.3 section 4.24, and it is **knots, then weights, then
        coordinates**::

            126, K, M, PROP1, PROP2, PROP3, PROP4,
            T(-M) .. T(K+1),          # K + M + 2 knots
            W(0) .. W(K),             # K + 1 weights (absent when polynomial)
            X(0), Y(0), Z(0), .. X(K), Y(K), Z(K)

        Reading the tail as ``control .. knots`` (which this method did until the
        KVLCC2/JBC audit in ``docs/baseline_fleet_type126_falsified.md``) put
        coordinates into the knot vector: on JBC DE 8 the "knots" came out as
        ``106.995, -3365.94, -0.0, …`` -- a non-monotone sequence of millimetres --
        and the "control points" as ``0.462, 0.462, 0.538, 0.538``, i.e. the knot
        vector itself.  The geometry never left the file; the reader indexed it wrong.
        """
        values = list(params)
        if int(values[0]) != ENTITY_BSPLINE_CURVE:
            raise IgesError(f"not a type-126 curve: {values[0]}")
        upper = int(values[1])  # K
        degree = int(values[2])  # M
        rational = int(values[6]) == 0  # PROP3: 0 = rational, 1 = polynomial
        count = upper + 1
        knot_count = count + degree + 1
        cursor = 7
        knots = tuple(values[cursor:cursor + knot_count])
        if len(knots) != knot_count:
            raise IgesError("type-126 curve truncated in its knot vector")
        cursor += knot_count
        if rational:
            weights = tuple(values[cursor:cursor + count])
            if len(weights) != count:
                raise IgesError("type-126 curve truncated in its weights")
            cursor += count
        else:
            weights = tuple(1.0 for _ in range(count))
        width = 3
        control = []
        for _ in range(count):
            chunk = values[cursor:cursor + width]
            cursor += width
            if len(chunk) != width:
                raise IgesError("type-126 curve truncated in its control points")
            control.append(tuple(chunk))
        return cls(
            degree=degree,
            knots=knots,
            control=tuple(control),
            weights=weights,
            dim=3,
        )


@dataclass(frozen=True)
class BSplineSurface:
    """Rational B-spline surface with independent (u, v) degrees (type 128)."""

    degree_u: int
    degree_v: int
    knots_u: tuple[float, ...]
    knots_v: tuple[float, ...]
    control: tuple[tuple[float, float, float], ...]  # row-major, u fastest
    weights: tuple[float, ...]
    count_u: int
    count_v: int
    #: Which control-point layout the reader chose ("interleaved" / "weight_block" /
    #: "polynomial"), and whether the choice was a coin flip between two plausible
    #: layouts.  A caller that integrates the surface (a hull mesher) can quarantine an
    #: ambiguous patch instead of silently meshing a mis-read control net.
    layout: str = "unknown"
    layout_ambiguous: bool = False

    @property
    def domain_u(self) -> tuple[float, float]:
        return self.knots_u[self.degree_u], self.knots_u[self.count_u]

    @property
    def domain_v(self) -> tuple[float, float]:
        return self.knots_v[self.degree_v], self.knots_v[self.count_v]

    def evaluate(self, u: float, v: float) -> tuple[float, float, float]:
        u_lo, u_hi = self.domain_u
        v_lo, v_hi = self.domain_v
        u = min(max(u, u_lo), u_hi)
        v = min(max(v, v_lo), v_hi)
        basis_u, span_u = _basis(self.count_u - 1, self.degree_u, self.knots_u, u)
        basis_v, span_v = _basis(self.count_v - 1, self.degree_v, self.knots_v, v)
        point = [0.0, 0.0, 0.0]
        denominator = 0.0
        for j in range(self.degree_v + 1):
            v_index = span_v - self.degree_v + j
            wv = basis_v[j]
            for i in range(self.degree_u + 1):
                u_index = span_u - self.degree_u + i
                flat = v_index * self.count_u + u_index
                weight = self.weights[flat] * basis_u[i] * wv
                xyz = self.control[flat]
                denominator += weight
                point[0] += weight * xyz[0]
                point[1] += weight * xyz[1]
                point[2] += weight * xyz[2]
        if denominator == 0.0:
            return (point[0], point[1], point[2])
        return (point[0] / denominator, point[1] / denominator, point[2] / denominator)

    def sample(self, count_u: int, count_v: int) -> list[tuple[float, float, float]]:
        """Row-major (v outer, u inner) samples on the *full* parametric domain."""
        if count_u < 2 or count_v < 2:
            raise ValueError("need at least two samples per direction")
        u_lo, u_hi = self.domain_u
        v_lo, v_hi = self.domain_v
        points = []
        for j in range(count_v):
            v = v_lo + (v_hi - v_lo) * j / (count_v - 1)
            for i in range(count_u):
                u = u_lo + (u_hi - u_lo) * i / (count_u - 1)
                points.append(self.evaluate(u, v))
        return points

    def bounding_box(self) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
        xs, ys, zs = zip(*self.control)
        return (min(xs), min(ys), min(zs)), (max(xs), max(ys), max(zs))

    @classmethod
    def from_parameters(cls, params: Sequence[float]) -> "BSplineSurface":
        """Build from a type-128 parameter list.

        Two control-point layouts occur in the wild and both are handled:

        * ``x, y, z, w`` interleaved per control point (IGES 5.3 as specified);
        * a block of ``n`` weights followed by a block of ``3n`` coordinates,
          which is what the NMRI/Tokyo-2015 writer emits.

        The layout is chosen by scoring the control net: the true net is smooth
        in both parametric directions, a mis-read one is not.
        """
        values = list(params)
        if int(values[0]) != ENTITY_BSPLINE_SURFACE:
            raise IgesError(f"not a type-128 surface: {values[0]}")
        upper_u, upper_v = int(values[1]), int(values[2])
        degree_u, degree_v = int(values[3]), int(values[4])
        rational = int(values[7]) == 0  # PROP3: 0 = rational, 1 = polynomial
        count_u, count_v = upper_u + 1, upper_v + 1
        points = count_u * count_v
        cursor = 10
        knot_u_count = count_u + degree_u + 1
        knot_v_count = count_v + degree_v + 1
        knots_u = tuple(values[cursor:cursor + knot_u_count])
        cursor += knot_u_count
        knots_v = tuple(values[cursor:cursor + knot_v_count])
        cursor += knot_v_count
        if len(knots_u) != knot_u_count or len(knots_v) != knot_v_count:
            raise IgesError("type-128 surface truncated in its knot vectors")
        body = values[cursor:]
        control, weights, layout, ambiguous = _split_control(
            body, points, rational, count_u, count_v
        )
        return cls(
            degree_u=degree_u,
            degree_v=degree_v,
            knots_u=knots_u,
            knots_v=knots_v,
            control=tuple(control),
            weights=tuple(weights),
            count_u=count_u,
            count_v=count_v,
            layout=layout,
            layout_ambiguous=ambiguous,
        )


def _net_roughness(
    control: Sequence[tuple[float, float, float]], count_u: int, count_v: int
) -> float:
    """Sum of squared second differences over the control net (lower = smoother)."""
    total = 0.0
    for j in range(count_v):
        base = j * count_u
        for i in range(1, count_u - 1):
            a, b, c = control[base + i - 1], control[base + i], control[base + i + 1]
            total += sum((a[k] - 2.0 * b[k] + c[k]) ** 2 for k in range(3))
    for i in range(count_u):
        for j in range(1, count_v - 1):
            a = control[(j - 1) * count_u + i]
            b = control[j * count_u + i]
            c = control[(j + 1) * count_u + i]
            total += sum((a[k] - 2.0 * b[k] + c[k]) ** 2 for k in range(3))
    return total


def _uniform_weights(weights: Sequence[float], tol: float = 1e-9) -> bool:
    """True when every weight is positive and equal to the others.

    Every writer seen in the workshop archives emits *unit* weights for the
    rational surfaces it writes (KCS, KVLCC2, JBC, ONRT all carry
    ``weights == 1.0`` throughout), so "this block is constant and positive" is a
    far stronger signal than "this net looks smooth" -- and it is the only signal
    that survives on a degenerate net, see :func:`_split_control`.
    """
    if not weights:
        return False
    lo, hi = min(weights), max(weights)
    if lo <= 0.0:
        return False
    return (hi - lo) <= tol * max(1.0, hi)


def _split_control(
    body: Sequence[float], points: int, rational: bool, count_u: int, count_v: int
) -> tuple[list[tuple[float, float, float]], list[float], str, bool]:
    """Split the control-point tail into points + weights, choosing the layout.

    Three layouts occur: IGES 5.3 interleaved ``x, y, z, w``; a block of weights
    followed by a block of coordinates (NMRI / Tokyo-2015 writer, which also
    appends the parameter bounds after the knots); and plain ``x, y, z`` for
    polynomial surfaces.  Writers also append a few trailing values that are not
    part of the control data, so the tail is matched from the front.

    **Layout choice.**  The original rule was "pick the net with the smallest sum
    of squared second differences".  That rule is undefined on a net with no
    interior control points: for a 2x2 or 2xn patch every second difference is
    identically zero, both candidate readings score 0.0, and the sort order (i.e.
    the string ``"interleaved"`` < ``"weight_block"``) silently decided.  Measured
    on the authoritative sources: 15 of the 132 ONRT patches and 2 of the 45 JBC
    patches were read this way, yielding control points inside the unit cube
    (``max|y| = 1.0 m`` on a hull whose half-beam is 9.39 m).  So the weight test
    comes first: a candidate whose weights are a *positive constant block* is
    preferred, because that is exactly what these writers emit.  The smoothness
    score is kept for genuinely rational patches with non-unit weights, where the
    weight test is silent.

    ``rational`` is the PROP3 flag (``False`` means rational, matching the IGES
    convention ``PROP3 = 0 ; 1 = polynomial``).  The polynomial layout is only a
    candidate when the record cannot be read as ``x, y, z, w`` quadruples or when
    PROP3 *declares* polynomial data: offering it unconditionally let a rational
    patch be parsed as truncated triples whenever the tail happened to be smooth.

    Returns ``(control, weights, layout, ambiguous)``.  ``ambiguous`` means the
    choice was not backed by the weight test and two layouts scored within a factor
    of four, so the caller (:func:`tensorcad.fleet_baseline.partition_surfaces`)
    may quarantine the patch instead of meshing a possibly mis-read net.
    """
    needed = 4 * points
    candidates: list[tuple[str, list[tuple[float, float, float]], list[float]]] = []
    if len(body) >= needed:
        interleaved = []
        weights_interleaved = []
        for index in range(points):
            x, y, z, w = body[4 * index:4 * index + 4]
            interleaved.append((x, y, z))
            weights_interleaved.append(w)
        candidates.append(("interleaved", interleaved, weights_interleaved))
        block = [tuple(body[points + 3 * index:points + 3 * index + 3]) for index in range(points)]
        candidates.append(("weight_block", block, list(body[:points])))  # type: ignore[arg-type]
    if len(body) >= 3 * points and (not rational or len(body) < needed):
        plain = [tuple(body[3 * index:3 * index + 3]) for index in range(points)]
        candidates.append(("polynomial", plain, [1.0] * points))  # type: ignore[arg-type]
    if not candidates:
        raise IgesError(
            f"control-point tail has {len(body)} values, too short for {points} points"
        )
    scored = []
    for name, control, weights in candidates:
        # A non-positive weight is not representable in these archives: every
        # rational patch here carries strictly positive weights, so a zero or
        # negative member means the tail was split at the wrong offset and a
        # coordinate was read as a weight.  Rejecting the candidate is what decides
        # a degenerate (2x2, 2xn) net, where the smoothness score is identically 0
        # for every reading -- that is how 15 ONRT patches kept a control point at
        # (1.0, 1.0, 1.0) m and spanned the unit cube.
        if any(not (w > 0.0) for w in weights) or any(len(point) != 3 for point in control):
            continue
        scored.append((_net_roughness(control, count_u, count_v), name, control, weights))
    if not scored:
        raise IgesError("no usable control-point layout for type-128 surface")
    uniform = [item for item in scored if _uniform_weights(item[3])]
    if len(uniform) == 1:
        best = uniform[0]
        return list(best[2]), list(best[3]), best[1], False
    if uniform:
        # Several readings claim constant weights.  ``weight_block`` read them out
        # of the record (``body[:points]``); ``polynomial`` merely assumes an
        # all-ones vector and swallows whatever constant block the record carries
        # as the first coordinates.  The block that is *in the data* wins: on a
        # 2x2 patch the polynomial reading is how the weight block ``1,1,1,1``
        # turned into the control point (1.0, 1.0, 1.0).
        from_record = [item for item in uniform if item[1] == "weight_block"]
        if len(from_record) == 1:
            best = from_record[0]
            return list(best[2]), list(best[3]), best[1], False
    scored.sort(key=lambda item: item[0])
    best = scored[0]
    ambiguous = False
    if len(scored) > 1 and best[0] > 0.0:
        second = scored[1][0]
        ambiguous = second > 0.0 and best[0] * 4.0 > second
    elif len(scored) > 1:
        # every candidate is perfectly smooth (degenerate net, no interior control
        # point): the smoothness test carries no information at all.
        ambiguous = True
    return list(best[2]), list(best[3]), best[1], ambiguous


# --------------------------------------------------------------------------- #
# directory entries and the file itself
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class DirectoryEntry:
    """One 20-field IGES directory entry."""

    sequence: int  # DE sequence number used by pointers
    entity_type: int
    parameter_pointer: int
    parameter_lines: int
    form: int
    label: str
    level: int
    status: str

    @property
    def blanked(self) -> bool:
        """Status field sub 1: blanked entities are not part of the model."""
        return self.status[0:1] == "0"


@dataclass
class Boundary:
    """Type-141 boundary: alternating model-space and parameter-space curves.

    Field order per IGES 5.3 (type 141): ``141, N, SPTR, N0, CTYPE(1..N)``.  The
    ``N`` curve pointers therefore start at parameter 5 (index 4); parameters 2
    and 3 carry the surface pointer and the parameter-space direction flag.  They
    are kept here because a caller that considers orientation needs them.
    """

    curves: tuple[int, ...]  # DE pointers, alternating (3D, (u,v), 3D, (u,v), ...)
    surface: int = 0
    parameter_flip: bool = False

    @classmethod
    def from_parameters(cls, params: Sequence[float]) -> "Boundary":
        count = int(params[1]) if len(params) > 1 else 0
        surface = int(params[2]) if len(params) > 2 else 0
        flip = bool(int(params[3])) if len(params) > 3 else False
        pointers = [int(p) for p in params[4:4 + count]]
        return cls(curves=tuple(pointers), surface=surface, parameter_flip=flip)


@dataclass
class CurveOnSurface:
    """Type-142 curve on a parametric surface: model-space curve + (u, v) loop.

    Field order per IGES 5.3 section 4.25 (type 142)::

        142, CRV, SPTR, BPTR, CPTR, IT

    ``SPTR`` is the basis surface, ``BPTR`` is the curve *in the parameter space of
    that surface* (the trim loop), ``CRV`` is the model-space curve, ``CPTR`` names a
    composite carrying the preferred representation of both, and ``IT`` selects that
    representation.

    This record was previously read one field to the right (``curve_3d`` = BPTR,
    ``curve_uv`` = CPTR), which made every trim loop a *model-space* polyline.  The
    measured decider is in ``IgesFile.trim_pointer_audit``: on the KVLCC2 archive the
    BPTR composite lies inside the basis surface's knot domain on 17 of 19 patches
    while the CPTR composite lies inside on none, so the loop is BPTR -- see
    ``artifacts/baseline_fleet/t126_trim_audit.json``.
    """

    surface: int
    curve_3d: int
    curve_uv: int | None
    preferred: int | None = None
    preferred_flag: int = 0

    @classmethod
    def from_parameters(cls, params: Sequence[float]) -> "CurveOnSurface":
        surface = int(params[2]) if len(params) > 2 else 0
        curve_3d = int(params[1]) if len(params) > 1 else 0
        curve_uv = int(params[3]) if len(params) > 3 else 0
        preferred = int(params[4]) if len(params) > 4 else 0
        flag = int(params[5]) if len(params) > 5 else 0
        return cls(
            surface=surface,
            curve_3d=curve_3d,
            curve_uv=curve_uv or None,
            preferred=preferred or None,
            preferred_flag=flag,
        )


@dataclass
class BoundedSurface:
    """Type-143/144 wrapper: basis surface plus its trimming boundaries.

    Field order follows IGES 5.3 section 4.26 (type 144): ``144, PTS, N1, N2,
    PTO, PTI``.  ``PTO`` -- the outer boundary -- is parameter **5** (index 4),
    and ``PTI`` -- the inner boundary -- is parameter 6 (index 5).  Reading the
    boundary out of index 3 (which carries ``N2``, normally ``0``) silently drops
    every trim loop, so a trimmed patch is then treated as its whole basis
    surface: that is a reader defect, not a property of the hull.
    """

    entity_type: int
    surface: int
    boundary: int | None
    inner: tuple[int, ...] = ()
    boundary_count: int = 0

    @classmethod
    def from_parameters(cls, params: Sequence[float]) -> "BoundedSurface":
        entity_type = int(params[0])
        surface = int(params[1])
        count = int(params[2]) if len(params) > 2 else 0
        outer = int(params[4]) if len(params) > 4 else 0
        inner = int(params[5]) if len(params) > 5 else 0
        return cls(
            entity_type=entity_type,
            surface=surface,
            boundary=outer or None,
            inner=(inner,) if inner else (),
            boundary_count=count,
        )


@dataclass
class IgesFile:
    """Parsed IGES file: global section plus entities Addressable by DE sequence."""

    path: str
    global_fields: tuple[str, ...]
    entries: dict[int, DirectoryEntry]
    raw_parameters: dict[int, list[float]]  # DE sequence -> numeric parameter list
    units: str
    unit_flag: int
    #: Raw units name (global field 15, e.g. ``"MM"``) and how the flag was decided
    #: (``"positional"`` / ``"name"`` / ``"name_overrides_flag"`` / ``"content_scan"``).
    unit_name: str = ""
    unit_source: str = "positional"
    #: Model-space scale (global field 13): multiply model coordinates by this to get
    #: real-space lengths.  IGES writes it in front of the units flag.
    model_scale: float = 1.0
    unit_warnings: tuple[str, ...] = ()
    #: Metres per model unit, model-space scale included.
    @property
    def metres_per_unit(self) -> float:
        return unit_scale(self.unit_flag) * self.model_scale

    # ---- lookups -------------------------------------------------------- #
    def of_type(self, entity_type: int) -> list[DirectoryEntry]:
        return [
            entry
            for _, entry in sorted(self.entries.items())
            if entry.entity_type == entity_type
        ]

    def parameters(self, sequence: int) -> list[float]:
        try:
            return self.raw_parameters[sequence]
        except KeyError as exc:
            raise IgesError(f"no parameter data for entity DE {sequence}") from exc

    def surfaces(self) -> list[BSplineSurface]:
        return [
            BSplineSurface.from_parameters(self.parameters(entry.sequence))
            for entry in self.of_type(ENTITY_BSPLINE_SURFACE)
        ]

    def surface_census(self) -> dict[str, object]:
        """What geometry the file carries, and how much of it the reader can bind.

        Entity identities are facts about IGES 5.3, not reader preferences:

        * **128** rational B-spline *surface* -- decoded into :class:`BSplineSurface`;
        * **126** rational B-spline *curve* -- decoded into :class:`BSplineCurve`;
          in the hull archives these are the trim boundaries (JBC 172, KVLCC2 304,
          ONRT 820), not dropped surfaces;
        * **141 / 143 / 144** boundary and trimmed-bounded patch wrappers;
        * **142** curve-on-surface, the parameter-space trim loop.

        Keys:
            ``surfaces_type128`` -- decodable basis surfaces;
            ``curves_type126`` -- rational B-spline curves;
            ``trimmed_patches_type143_144`` -- wrapper entities present;
            ``boundaries_type141`` / ``curves_on_surface_type142`` -- trim carriers;
            ``basis_surfaces_trimmed`` / ``basis_surfaces_untrimmed`` -- basis
            surfaces a wrapper names / that no wrapper names;
            ``by_type`` -- the full histogram.

        The last pair is the one a volume integrator needs: a hull whose basis
        surfaces are mostly *untrimmed* integrates over extra sheet, and that -- not
        a mis-decoded entity type -- is what a wrong displacement describes.
        """
        by_type: dict[int, float] = {}
        for entry in self.entries.values():
            by_type[entry.entity_type] = by_type.get(entry.entity_type, 0) + 1
        basis_total = by_type.get(ENTITY_BSPLINE_SURFACE, 0)
        trimmed: set[int] = set()
        for kind in (ENTITY_TRIMMED_SURFACE, ENTITY_BOUNDED_SURFACE):
            for entry in self.of_type(kind):
                wrapper = BoundedSurface.from_parameters(self.parameters(entry.sequence))
                target = self.resolve(wrapper.surface) if wrapper.surface else None
                if target is not None and target.entity_type == ENTITY_BSPLINE_SURFACE:
                    trimmed.add(target.sequence)
        return {
            "surfaces_type128": basis_total,
            "curves_type126": by_type.get(ENTITY_BSPLINE_CURVE, 0),
            "trimmed_patches_type143_144": by_type.get(ENTITY_BOUNDED_SURFACE, 0)
            + by_type.get(ENTITY_TRIMMED_SURFACE, 0),
            "boundaries_type141": by_type.get(ENTITY_BOUNDARY, 0),
            "curves_on_surface_type142": by_type.get(ENTITY_CURVE_ON_SURFACE, 0),
            "basis_surfaces_trimmed": len(trimmed),
            "basis_surfaces_untrimmed": basis_total - len(trimmed),
            "by_type": dict(sorted(by_type.items())),
        }

    def trim_pointer_audit(self, samples: int = 40) -> dict[str, object]:
        """Decide, per trimmed patch, which curve-on-surface pointer lives in (u, v).

        Type 142 carries two curve pointers (the model-space curve and the
        parameter-space loop) plus the preferred-representation composite, and the
        *field order alone has not been enough to tell them apart on real bytes*:
        the shape of the record is identical either way, and the writer's directory
        numbering (first vs second record) is itself ambiguous.  So the question is
        answered by measurement instead of by assertion: a parameter-space loop must
        lie inside the basis surface's knot domain, a model-space curve does not.

        For every trimmed patch this samples both candidate composites, reports the
        fraction of samples inside ``domain_u x domain_v``, and marks a patch
        ``decided`` only when one candidate is inside on at least half its samples
        and leads the other by more than 10 points.  Patches that are not decided are
        reported, never defaulted -- a caller that meshes a trim loop has to know.

        Keys: ``patches`` (row per wrapper), ``decided_field_crv`` (rows where the
        ``curve_3d`` field won), ``decided_field_bptr`` (rows where ``curve_uv`` won),
        ``undecided`` (rows neither candidate claims), ``samples_per_candidate``.
        """
        rows: list[dict[str, object]] = []
        decided_3d = decided_uv = undecided = 0
        for index, (wrapper, basis, curves) in enumerate(self.trimmed_surfaces()):
            if basis is None or not curves:
                rows.append(
                    {
                        "patch": index,
                        "surface_de": wrapper.surface,
                        "basis_bound": basis is not None,
                        "curves": len(curves),
                        "verdict": "no_basis_or_no_curve",
                    }
                )
                undecided += 1
                continue
            u0, u1 = basis.domain_u
            v0, v1 = basis.domain_v
            fractions: dict[str, float | None] = {}
            for label, pointer in (("field_crv", curves[0].curve_3d), ("field_bptr", curves[0].curve_uv)):
                fractions[label] = self._inside_fraction(pointer, u0, u1, v0, v1, samples)
            best = None
            for label, value in fractions.items():
                if value is None or value < 0.5:
                    continue
                if best is None or value > (fractions[best] or 0.0):
                    best = label
            margin = None
            if best is not None:
                other = "field_bptr" if best == "field_crv" else "field_crv"
                margin = round((fractions[best] or 0.0) - (fractions[other] or 0.0), 3)
                if margin <= 0.1:
                    best, margin = None, margin
            verdict = "decided" if best else "undecided"
            if best == "field_crv":
                decided_3d += 1
            elif best == "field_bptr":
                decided_uv += 1
            else:
                undecided += 1
            rows.append(
                {
                    "patch": index,
                    "surface_de": wrapper.surface,
                    "basis_bound": True,
                    "curves": len(curves),
                    "domain": [[round(u0, 6), round(u1, 6)], [round(v0, 6), round(v1, 6)]],
                    "inside_fraction": {
                        key: (None if value is None else round(value, 3)) for key, value in fractions.items()
                    },
                    "winner_field": best,
                    "margin": margin,
                    "verdict": verdict,
                }
            )
        return {
            "patches": rows,
            "patch_count": len(rows),
            "decided_field_crv": decided_3d,
            "decided_field_bptr": decided_uv,
            "undecided": undecided,
            "samples_per_candidate": samples,
        }

    def _inside_fraction(
        self,
        pointer: int | None,
        u0: float,
        u1: float,
        v0: float,
        v1: float,
        samples: int = 40,
    ) -> float | None:
        """Fraction of a composite's samples that fall inside the (u, v) domain."""
        if not pointer:
            return None
        points: list[tuple[float, ...]] = []
        try:
            entity = self.resolve(pointer)
            if entity is None:
                return None
            if entity.entity_type == ENTITY_COMPOSITE_CURVE:
                params = self.parameters(entity.sequence)
                count = int(params[1])
                for child in params[2:2 + count]:
                    target = self.resolve(int(child))
                    if target is None or target.entity_type != ENTITY_BSPLINE_CURVE:
                        continue
                    curve = BSplineCurve.from_parameters(self.parameters(target.sequence))
                    points.extend(curve.sample(samples))
            elif entity.entity_type == ENTITY_BSPLINE_CURVE:
                curve = BSplineCurve.from_parameters(self.parameters(entity.sequence))
                points.extend(curve.sample(samples))
        except (IgesError, ValueError, IndexError):
            return None
        if not points:
            return None
        inside = 0
        for point in points:
            u, v = point[0], point[1]
            if u0 - 1e-9 <= u <= u1 + 1e-9 and v0 - 1e-9 <= v <= v1 + 1e-9:
                inside += 1
        return inside / len(points)

    def resolve(self, pointer: int) -> DirectoryEntry | None:
        """Look up the entity a pointer carried by another entity refers to.

        Writers disagree about which of a directory entry's two 80-column records
        the pointer names.  IGES 5.3 numbers the *first* record, and the workshop
        archives split both ways: the ONRT/Rhino export carries odd pointers
        (``141,1,3,27,...`` names DE 27) while this reader keys every entry by its
        second record (``second[73:80]``), i.e. the same entry is key 28.  A lookup
        by exact key therefore missed every reference in that file, and
        ``trimmed_surfaces`` failed with "boundary DE 1 missing" instead of
        resolving the patch.

        The neighbour is only consulted for an *odd* pointer, whose second record
        is ``pointer + 1`` -- an entry cannot have a third record, so the step can
        never land on an unrelated entity.
        """
        entry = self.entries.get(pointer)
        if entry is not None:
            return entry
        if pointer % 2 == 1:
            return self.entries.get(pointer + 1)
        return None

    def curve(self, sequence: int) -> BSplineCurve:
        params = self.parameters(sequence)
        entity_type = int(params[0])
        if entity_type != ENTITY_BSPLINE_CURVE:
            raise IgesError(f"DE {sequence} is type {entity_type}, not a B-spline curve")
        return BSplineCurve.from_parameters(params)

    def trim_polygon_uv(self, pointer: int, samples: int = 200) -> list[tuple[float, float]]:
        """Sample a trim boundary into a closed polygon in the surface's (u, v) space.

        ``pointer`` is a DE pointer to a type-126 curve (the loop directly) or to a
        type-102 composite curve whose children are type-126 curves (a loop written
        as several pieces).  Model-space curves are *not* valid here: a trim loop
        has to be tested in the parameter space of the surface it trims.
        """
        entry = self.resolve(pointer)
        if entry is None:
            raise IgesError(f"trim boundary DE {pointer} missing")
        if entry.entity_type == ENTITY_COMPOSITE_CURVE:
            params = self.parameters(entry.sequence)
            count = int(params[1])
            polygon: list[tuple[float, float]] = []
            for child in params[2:2 + count]:
                polygon.extend(self.trim_polygon_uv(int(child), samples=max(8, samples // max(count, 1))))
            return polygon
        if entry.entity_type != ENTITY_BSPLINE_CURVE:
            raise IgesError(
                f"trim boundary DE {pointer} is type {entry.entity_type}, not a B-spline/composite curve"
            )
        curve = BSplineCurve.from_parameters(self.parameters(entry.sequence))
        return [(point[0], point[1]) for point in curve.sample(samples)]

    def trimmed_surfaces(self) -> list[tuple[BoundedSurface, BSplineSurface | None, list["CurveOnSurface"]]]:
        """Return ``(wrapper, basis surface, boundary curves)`` for each trimmed patch."""
        result = []
        for entry in self.of_type(ENTITY_TRIMMED_SURFACE) + self.of_type(ENTITY_BOUNDED_SURFACE):
            wrapper = BoundedSurface.from_parameters(self.parameters(entry.sequence))
            basis = None
            basis_entry = self.resolve(wrapper.surface)
            if basis_entry is not None and basis_entry.entity_type == ENTITY_BSPLINE_SURFACE:
                basis = BSplineSurface.from_parameters(self.parameters(basis_entry.sequence))
            curves: list[CurveOnSurface] = []
            if wrapper.boundary:
                boundary_entry = self.resolve(wrapper.boundary)
                if boundary_entry is None:
                    raise IgesError(f"boundary DE {wrapper.boundary} missing")
                if boundary_entry.entity_type == ENTITY_BOUNDARY:
                    boundary = Boundary.from_parameters(self.parameters(boundary_entry.sequence))
                    for pointer in boundary.curves:
                        target = self.resolve(pointer)
                        if target is not None and target.entity_type == ENTITY_CURVE_ON_SURFACE:
                            curves.append(CurveOnSurface.from_parameters(self.parameters(target.sequence)))
                elif boundary_entry.entity_type == ENTITY_CURVE_ON_SURFACE:
                    curves.append(
                        CurveOnSurface.from_parameters(self.parameters(boundary_entry.sequence))
                    )
                for pointer in wrapper.inner:
                    target = self.resolve(pointer)
                    if target is not None and target.entity_type == ENTITY_CURVE_ON_SURFACE:
                        curves.append(CurveOnSurface.from_parameters(self.parameters(target.sequence)))
            result.append((wrapper, basis, curves))
        return result

    def stats(self) -> dict[str, int]:
        counts: dict[str, int] = {}
        for entry in self.entries.values():
            counts[str(entry.entity_type)] = counts.get(str(entry.entity_type), 0) + 1
        return dict(sorted(counts.items(), key=lambda item: int(item[0])))


def _split_sections(lines: Sequence[str]) -> dict[str, list[str]]:
    sections: dict[str, list[str]] = {letter: [] for letter in _SECTIONS}
    for line in lines:
        if len(line) < 73:
            continue
        letter = line[72]
        if letter in sections:
            sections[letter].append(line)
    if not sections["D"]:
        raise IgesError("file has no directory section; not an IGES 5.x file")
    return sections


#: Metres per unit for the IGES global unit flag.
UNIT_SCALE: dict[int, float] = {
    1: 0.0254,  # inches
    2: 0.001,  # millimetres
    3: 0.01,  # centimetres
    4: 0.0254,  # inches
    5: 0.3048,  # feet
    6: 1.0,  # metres
    7: 1000.0,  # kilometres
    8: 2.54e-5,  # mils
    9: 1e-6,  # microns
    10: 0.1,  # decimetres
    11: 100.0,  # hectometres
}


def unit_scale(unit_flag: int) -> float:
    """Metres per model unit for an IGES global unit flag (unknown -> metres)."""
    return UNIT_SCALE.get(unit_flag, 1.0)


#: IGES units-name (global field 15) -> units flag (global field 14), for the flags
#: :data:`UNIT_SCALE` knows how to convert.  The name is the field a human reads off the
#: file, so it is used to cross-check (and, on disagreement, to correct) the numeric flag.
UNIT_NAME_FLAGS: dict[str, int] = {
    "IN": 1, "MM": 2, "CM": 3, "INCH": 1, "FT": 5, "M": 6, "KM": 7,
    "MIL": 8, "UM": 9, "DM": 10, "HM": 11,
}


def hollerith_fields(payload: str) -> tuple[str, ...]:
    """Split the IGES global section into its numbered fields.

    The global section is comma separated *except* that a non-empty field is a
    Hollerith string whose payload may itself contain commas, and the two delimiter
    fields at the front (``1H,`` and ``1H;``) are those delimiter characters written as
    Hollerith data.  A plain ``str.split(",")`` therefore yields a different number of
    tokens depending on whether a writer emits those two fields, which shifts every
    later field index.  That is exactly what used to mis-read the units flag: the KCS
    file leaves fields 1-2 empty (unit flag on token 13) while the Rhino files
    JBC/ONRT write ``1H,``/``1H;`` (unit flag on token 14), and the old content scan
    returned 6 ("metres") for JBC, whose file is millimetres -- a silent factor of 1000.

    The tokenizer is the tested one from :mod:`tensorcad.geomquality.iges_trim`; it is
    imported lazily so the dependency stays one-way.
    """
    from .geomquality.iges_trim import _tokenize_parameters

    return tuple(_tokenize_parameters(payload))


def _global_units(
    fields: Sequence[str],
) -> tuple[int, str, float, str, list[str]]:
    """Return ``(unit_flag, unit_name, model_space_scale, source, warnings)``.

    Positional first (global fields 13 model-space scale, 14 units flag, 15 units
    name), then corrected by the units *name* when the two disagree, and only then by
    the legacy content scan (which is unreliable: fields 7-11, the precision
    magnitudes, are a real number followed by a small integer, so they mimic the
    scale/flag pair).
    """
    warnings: list[str] = []
    flag: int | None = None
    name = ""
    model_scale = 1.0
    source = "positional"
    if len(fields) > 13:
        try:
            candidate = int(float(fields[13].strip()))
        except ValueError:
            candidate = 0
        if 1 <= candidate <= 11:
            flag = candidate
        if len(fields) > 14:
            name = fields[14].strip()
    named = UNIT_NAME_FLAGS.get(name.upper(), None) if name else None
    if named is not None:
        if flag is not None and flag != named:
            warnings.append(
                f"global units flag {flag} disagrees with units name {name!r}: "
                f"the name ({named}) was used"
            )
            source = "name_overrides_flag"
        elif flag is None:
            source = "name"
        flag = named
    if flag is None:
        # Last resort: the historical content scan, kept for files with a truncated
        # global section.  It is what mis-read JBC, so it is only used when nothing
        # else identifies the units, and it is flagged as such.
        flag = _unit_flag_content_scan(fields)
        source = "content_scan"
    try:
        model_scale = float(fields[12].strip().replace("D", "E").replace("d", "e"))
    except (ValueError, IndexError):
        model_scale = 1.0
    if not (model_scale > 0.0):
        warnings.append(f"model-space scale {fields[12]!r} is not positive: using 1.0")
        model_scale = 1.0
    return flag, name, model_scale, source, warnings


def _unit_flag_content_scan(fields: Sequence[str]) -> int:
    """Historical fallback: first real number followed by a small integer."""
    for index, field in enumerate(fields[:-1]):
        text = field.strip().replace("D", "E").replace("d", "e")
        try:
            float(text)
        except ValueError:
            continue
        following = fields[index + 1].strip()
        if following.isdigit() and 1 <= int(following) <= 11:
            return int(following)
    if len(fields) > 13 and fields[13].strip().isdigit():
        return int(fields[13].strip())
    return 1


def _unit_flag(fields: Sequence[str]) -> int:
    """Units flag for a tokenised global section (see :func:`_global_units`)."""
    return _global_units(fields)[0]


def read_iges(path: str) -> IgesFile:
    """Read ``path`` and return the parsed model."""
    with open(path, "r", errors="replace") as handle:
        lines = handle.read().splitlines()
    sections = _split_sections(lines)

    global_payload = "".join(line[0:72] for line in sections["G"])
    global_fields = hollerith_fields(global_payload)
    unit_flag, unit_name, model_scale, unit_source, unit_warnings = _global_units(global_fields)
    units = {1: "in", 2: "mm", 3: "cm", 4: "in", 5: "ft", 6: "m", 7: "km"}.get(unit_flag, "?")

    parameter_records: dict[int, str] = {}
    for line in sections["P"]:
        try:
            sequence = int(line[73:80])
        except ValueError:  # pragma: no cover - malformed input
            continue
        # Columns 1-64 carry the free-form data per IGES 5.3, but some writers
        # (the SIMMAN 2020 ONRT file) let a field run past column 64 and break the
        # number at column 72. When a record does not end on a separator the
        # overflow columns belong to the same token, so they are re-appended
        # instead of being dropped (which used to split "1.0D0" in half).
        data = line[0:64]
        overflow = line[64:72]
        if data.rstrip()[-1:] not in (",", ";") and overflow.strip():
            data += overflow
        parameter_records[sequence] = data

    entries: dict[int, DirectoryEntry] = {}
    for index in range(0, len(sections["D"]) - 1, 2):
        first, second = sections["D"][index], sections["D"][index + 1]
        try:
            entity_type = int(first[0:8])
            parameter_pointer = int(first[8:16])
            our_sequence = int(second[73:80])
            parameter_lines = int(second[24:32].strip() or 1)
            form = int(second[40:48].strip() or 0)
            level = int(first[32:40].strip() or 0)
        except ValueError as exc:
            raise IgesError(f"unreadable directory entry at line {index}") from exc
        entries[our_sequence] = DirectoryEntry(
            sequence=our_sequence,
            entity_type=entity_type,
            parameter_pointer=parameter_pointer,
            parameter_lines=parameter_lines,
            form=form,
            label=second[56:64].strip(),
            level=level,
            status=first[64:72].strip(),
        )

    raw_parameters: dict[int, list[float]] = {}
    for sequence, entry in entries.items():
        payload = ""
        for offset in range(max(entry.parameter_lines, 1)):
            record = parameter_records.get(entry.parameter_pointer + offset)
            if record is None:
                break
            payload += record
        if not payload:
            continue
        payload = payload.split(";")[0]
        values: list[float] = []
        for token in _tokens(payload):
            try:
                values.append(_strip_number(token))
            except IgesError:
                # Hollerith strings (names, ledger text on type 314/406 entities
                # and the workspace notes some writers attach) are not numbers.
                # They are replaced by NaN so the *positions* of the numeric
                # fields stay stable for every entity this module does parse.
                values.append(float("nan"))
        raw_parameters[sequence] = values

    return IgesFile(
        path=path,
        global_fields=global_fields,
        entries=entries,
        raw_parameters=raw_parameters,
        units=units,
        unit_flag=unit_flag,
        unit_name=unit_name,
        unit_source=unit_source,
        model_scale=model_scale,
        unit_warnings=tuple(unit_warnings),
    )


def read_iges_surfaces(path: str) -> list[BSplineSurface]:
    """Convenience wrapper: every type-128 surface in the file."""
    return read_iges(path).surfaces()