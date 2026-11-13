"""Compatible edge representations: equal degree **and** equal knot vector.

Degree elevation (:meth:`tensorcad.nurbs.NURBSCurve.elevate_degree` /
:meth:`tensorcad.nurbs.NURBSSurface.elevate_degree`) removes one of the two
mismatches between patches that share an edge. It does not remove the other: two
edges of the same degree can still disagree in their **knot vector**, and then they
carry different numbers of control points for the same curve, so a stitcher has no
shared data to weld -- it has to fall back to a numeric coincidence test inside some
tolerance. This module lands the missing half:

* :func:`insert_knot` -- Boehm's single-knot insertion (Piegl & Tiller A5.1) on a
  homogeneous control polygon, so it is exact for rational curves as well;
* :func:`align_edges` -- elevate both edges to the common degree, then insert knots
  until the two knot vectors are *equal*, which makes two edges that describe the
  same curve come out with identical control points and weights: a watertight seam
  by construction rather than by tolerance;
* :func:`boundary_curve` -- the four boundary curves of a patch, so the seam of two
  patches can be taken without any index arithmetic in the caller;
* :func:`edge_report` -- the measured state of a seam: degrees, control counts,
  knot-vector equality, control-point/weight difference and the geometric gap at the
  same parameters. It is what makes the claim checkable, and it also reports the gap
  of a seam that does **not** coincide: alignment makes representations compatible,
  it never makes geometry coincide.

Named limitations
-----------------
* Seams, not assemblies: :func:`align_patch_seam` refines the two patches that own a
  shared edge until their boundary rows are identical, so the seam is watertight by
  construction -- in the patch data, not only in a pair of loose curves. What stays
  out of scope is the *assembly*: no topology, no trimming, no orientation search over
  a whole network, and a patch that is a neighbour of two refined patches is refined
  twice (which is exact, just not minimal).
* The two edges must share a parameter domain. Insertion can equalise the *knot
  vector*, it cannot equalise the *parametrisation*: two edges of the same circle
  split into three and into six arcs have the same point set but visit it at different
  parameters, and after alignment they still disagree at the same ``u``
  (``examples/nurbs_edge_alignment`` measures 3.85e-2 m on a 1 m radius, 3.9% of the
  radius, with weights 0.75 against ``cos 30 deg``). That is a different curve in the
  sense a solver cares about, and the report says so instead of hiding it behind "the
  representations now match".
* Knot insertion stops at multiplicity ``degree`` (a ``C0`` break); a higher
  multiplicity would split the curve, so it is refused and the fix (elevate first) is
  named.
* Direction handling is explicit: :func:`reverse_curve` reparametrises an edge whose
  neighbour runs the other way; :func:`align_edges` does not guess it.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass

from .highorder.nurbs_geom import insert_knot as _insert_surface_knot
from .nurbs import (
    NURBSCurve,
    NURBSError,
    NURBSSurface,
    _dehomogenise,
    _rescale_to_uniform,
    _uniform_weight,
    find_span,
)

__all__ = [
    "SIDES",
    "SeamAlignment",
    "align_edges",
    "align_patch_seam",
    "boundary_curve",
    "edge_report",
    "insert_knot",
    "reverse_curve",
]

SIDES = ("v0", "v1", "u0", "u1")

_DOMAIN_TOLERANCE = 1e-12


# -- knot insertion ---------------------------------------------------------
def _homogeneous_points(curve: NURBSCurve) -> list[tuple[float, float, float, float]]:
    return [
        (weight * point[0], weight * point[1], weight * point[2], weight)
        for point, weight in zip(curve.control_points, curve.weights, strict=True)
    ]


def _multiplicity(knots: Sequence[float], value: float) -> int:
    return sum(1 for knot in knots if knot == value)


def insert_knot(curve: NURBSCurve, value: float, times: int = 1) -> NURBSCurve:
    """Insert the knot ``value`` ``times`` times into ``curve`` (Boehm A5.1).

    The geometry and the parametrisation are unchanged: the returned curve evaluates
    to the same point at the same parameter, and its domain is the same. What changes
    is the representation -- one more control point per insertion -- which is the
    degree of freedom that lets two edges be given the same knot vector.

    The insertion is done on the homogeneous control polygon, so a rational curve
    stays exactly the same rational curve (inserting into the 3-D net and carrying
    the weights across is a different, generally wrong curve).
    """
    if not isinstance(times, int) or isinstance(times, bool):
        raise NURBSError(f"times must be a non-negative integer, got {times!r}")
    if times < 0:
        raise NURBSError(f"times must be a non-negative integer, got {times!r}")
    if times == 0:
        return curve

    degree = curve.degree
    low, high = curve.domain
    if not (low <= value <= high):
        raise NURBSError(f"knot {value!r} is outside the domain [{low!r}, {high!r}]")
    existing = _multiplicity(curve.knots, value)
    if existing + times > degree:
        raise NURBSError(
            f"knot {value!r} already has multiplicity {existing}; inserting {times} "
            f"more would exceed degree {degree} and break the curve there -- elevate "
            "the degree first if the curve needs more control there"
        )

    points = _homogeneous_points(curve)
    knots = list(curve.knots)
    for _ in range(times):
        span = find_span(tuple(knots), degree, value)
        head = points[: span - degree + 1]
        blended: list[tuple[float, float, float, float]] = []
        for index in range(span - degree + 1, span + 1):
            denominator = knots[index + degree] - knots[index]
            alpha = 0.0 if denominator == 0.0 else (value - knots[index]) / denominator
            left = points[index]
            right = points[index - 1]
            blended.append(
                (
                    alpha * left[0] + (1.0 - alpha) * right[0],
                    alpha * left[1] + (1.0 - alpha) * right[1],
                    alpha * left[2] + (1.0 - alpha) * right[2],
                    alpha * left[3] + (1.0 - alpha) * right[3],
                )
            )
        points = head + blended + points[span:]
        knots.insert(span + 1, value)

    uniform = _uniform_weight(points)
    # Keep an all-equal weight net exactly all-equal: the blend returns ``w +/- 1e-16``,
    # which would silently move the curve onto the rational evaluation path with no
    # geometric change but with a changed *type*.
    rescaled = points if uniform is None else _rescale_to_uniform(points, uniform)
    weights: tuple[float, ...] = (
        tuple(point[3] for point in rescaled)
        if uniform is None
        else tuple(uniform for _ in rescaled)
    )
    return NURBSCurve(
        degree=degree,
        knots=tuple(knots),
        control_points=tuple(_dehomogenise(point) for point in rescaled),
        weights=weights,
    )


# -- edges ------------------------------------------------------------------
def boundary_curve(surface: NURBSSurface, side: str, *, reverse: bool = False) -> NURBSCurve:
    """One boundary curve of a patch: ``"v0"``, ``"v1"``, ``"u0"`` or ``"u1"``.

    ``"v0"`` is the edge at ``v = v_min`` (a curve in ``u``), ``"u0"`` the edge at
    ``u = u_min`` (a curve in ``v``). ``reverse=True`` returns the same edge traced
    the other way, which is what a neighbour patch that runs backwards needs before
    the two can share control points.
    """
    if side not in SIDES:
        raise ValueError(f"side must be one of {SIDES}, got {side!r}")
    degree_u, degree_v = surface.degree
    net = surface.control_net
    weights = surface.weights
    if side == "v0":
        curve = NURBSCurve(
            degree=degree_u,
            knots=surface.knots[0],
            control_points=tuple(row[0] for row in net),
            weights=tuple(row[0] for row in weights),
        )
    elif side == "v1":
        curve = NURBSCurve(
            degree=degree_u,
            knots=surface.knots[0],
            control_points=tuple(row[-1] for row in net),
            weights=tuple(row[-1] for row in weights),
        )
    elif side == "u0":
        curve = NURBSCurve(
            degree=degree_v,
            knots=surface.knots[1],
            control_points=tuple(net[0]),
            weights=tuple(weights[0]),
        )
    else:
        curve = NURBSCurve(
            degree=degree_v,
            knots=surface.knots[1],
            control_points=tuple(net[-1]),
            weights=tuple(weights[-1]),
        )
    return reverse_curve(curve) if reverse else curve


def reverse_curve(curve: NURBSCurve) -> NURBSCurve:
    """The same edge traced in the opposite direction (knots mirrored onto the domain).

    Reversing the control polygon alone would not reverse the parametrisation: the
    knot vector has to be reflected about the domain midpoint so that ``u`` of the
    result visits the original at ``low + high - u``.
    """
    low, high = curve.domain
    knots = tuple(low + high - knot for knot in reversed(curve.knots))
    return NURBSCurve(
        degree=curve.degree,
        knots=knots,
        control_points=tuple(reversed(curve.control_points)),
        weights=tuple(reversed(curve.weights)),
    )


# -- alignment --------------------------------------------------------------
def _interior_multiplicities(curve: NURBSCurve) -> dict[float, int]:
    low, high = curve.domain
    values: dict[float, int] = {}
    for knot in curve.knots:
        if low < knot < high:
            values[knot] = values.get(knot, 0) + 1
    return values


def align_edges(first: NURBSCurve, second: NURBSCurve) -> tuple[NURBSCurve, NURBSCurve]:
    """Give two edges one representation: common degree, common knot vector.

    Steps, each forced by the next:

    1. **Common degree** -- elevate the lower-degree edge by the difference
       (:meth:`NURBSCurve.elevate_degree`), which grows every multiplicity by the same
       amount and therefore keeps each edge's own continuity.
    2. **Common breakpoints** -- a knot value that one edge carries and the other
       does not is inserted into the other, so the two end up with the same set of
       breakpoints.
    3. **Common multiplicities** -- each shared breakpoint is raised to the larger of
       the two multiplicities (:func:`insert_knot`).

    Afterwards the two knot vectors are equal, so two edges that describe the same
    curve have the same number of control points and the same values: the seam is
    watertight by construction. Two edges that describe *different* curves keep a gap
    -- :func:`edge_report` measures it, and the alignment does not hide it.
    """
    low_a, high_a = first.domain
    low_b, high_b = second.domain
    if not (
        math.isclose(low_a, low_b, rel_tol=_DOMAIN_TOLERANCE, abs_tol=_DOMAIN_TOLERANCE)
        and math.isclose(high_a, high_b, rel_tol=_DOMAIN_TOLERANCE, abs_tol=_DOMAIN_TOLERANCE)
    ):
        raise NURBSError(
            "align_edges needs a common parameter domain; got "
            f"[{low_a!r}, {high_a!r}] and [{low_b!r}, {high_b!r}] -- reparametrise "
            "the edges onto one range first"
        )

    target = max(first.degree, second.degree)
    left = first.elevate_degree(target - first.degree)
    right = second.elevate_degree(target - second.degree)

    # Union of breakpoints: a knot value one edge carries and the other does not is
    # inserted into the other (that is exactly what insertion is for), and a value
    # both carry is raised to the larger multiplicity. After this the two knot
    # vectors are equal, so two edges describing the same parametrised curve have the
    # same control points and weights.
    left_knots = _interior_multiplicities(left)
    right_knots = _interior_multiplicities(right)
    for value in sorted(set(left_knots) | set(right_knots)):
        wanted = max(left_knots.get(value, 0), right_knots.get(value, 0))
        deficit = wanted - left_knots.get(value, 0)
        if deficit:
            left = insert_knot(left, value, deficit)
        deficit = wanted - right_knots.get(value, 0)
        if deficit:
            right = insert_knot(right, value, deficit)
    return left, right


# -- measurement ------------------------------------------------------------
def _sample_parameters(first: NURBSCurve, second: NURBSCurve, samples: int) -> tuple[float, ...]:
    low, high = first.domain
    values = {low + (high - low) * index / (samples - 1) for index in range(samples)}
    for curve in (first, second):
        for knot in curve.knots:
            if low <= knot <= high:
                values.add(knot)
    return tuple(sorted(values))


def _max_distance(
    left: tuple[float, float, float], right: tuple[float, float, float]
) -> float:
    return math.dist(left, right)


def edge_report(
    first: NURBSCurve,
    second: NURBSCurve,
    *,
    samples: int = 129,
    tolerance: float = 1e-12,
) -> dict[str, object]:
    """Measure a seam: representation agreement, and the geometry it is meant to carry.

    ``max_gap`` is the geometric distance at the **same parameters**, which is the only
    seam measure that can distinguish "the representation now matches" from "the
    geometry now matches". ``watertight`` combines both: equal degree, equal knot
    vector, and control points plus weights agreeing to ``tolerance``.
    """
    parameters = _sample_parameters(first, second, samples)
    gap = max(_max_distance(first.evaluate(u), second.evaluate(u)) for u in parameters)

    same_degree = first.degree == second.degree
    knots_equal = same_degree and first.knots == second.knots
    same_count = len(first.control_points) == len(second.control_points)
    control_difference = (
        max(
            _max_distance(point_a, point_b)
            for point_a, point_b in zip(
                first.control_points, second.control_points, strict=True
            )
        )
        if same_count
        else None
    )
    weight_difference = (
        max(
            abs(weight_a - weight_b)
            for weight_a, weight_b in zip(first.weights, second.weights, strict=True)
        )
        if same_count
        else None
    )
    return {
        "degree": [first.degree, second.degree],
        "control_points": [len(first.control_points), len(second.control_points)],
        "rational": [first.rational, second.rational],
        "domain": [list(first.domain), list(second.domain)],
        "same_degree": same_degree,
        "knots_equal": knots_equal,
        "max_control_difference": control_difference,
        "max_weight_difference": weight_difference,
        "max_gap": gap,
        "tolerance": tolerance,
        "watertight": bool(
            knots_equal
            and control_difference is not None
            and control_difference <= tolerance
            and weight_difference is not None
            and weight_difference <= tolerance
        ),
    }


# -- patch-level seams ------------------------------------------------------
_EDGE_AXIS = {"v0": 0, "v1": 0, "u0": 1, "u1": 1}


@dataclass(frozen=True)
class SeamAlignment:
    """Two refined patches plus the evidence that the seam is (or is not) closed.

    ``first``/``second`` are the patches to use from now on -- carrying the same
    geometry as the inputs, with their shared boundary rows made identical.
    ``report`` is :func:`edge_report` on the two boundary curves after the refinement,
    extended with the refinement's own bookkeeping, and ``operations`` is the audit
    trail of what was done to which patch (degree elevation and knot insertions).
    """

    first: NURBSSurface
    second: NURBSSurface
    first_edge: NURBSCurve
    second_edge: NURBSCurve
    report: dict[str, object]
    operations: tuple[dict[str, object], ...]

    @property
    def watertight(self) -> bool:
        """Whether the two patches now share one edge representation."""
        return bool(self.report["watertight"])

    def to_dict(self) -> dict[str, object]:
        return {
            "report": dict(self.report),
            "operations": [dict(entry) for entry in self.operations],
            "first_edge": {
                "degree": self.first_edge.degree,
                "control_points": len(self.first_edge.control_points),
            },
            "second_edge": {
                "degree": self.second_edge.degree,
                "control_points": len(self.second_edge.control_points),
            },
        }


def _axis_multiplicities(knots: Sequence[float], degree: int, count: int) -> dict[float, int]:
    """Interior multiplicities of a clamped knot vector with ``count`` basis functions."""
    low, high = knots[degree], knots[count]
    values: dict[float, int] = {}
    for knot in knots:
        if low < knot < high:
            values[knot] = values.get(knot, 0) + 1
    return values


def _refine_axis(
    surface: NURBSSurface,
    axis: int,
    target: NURBSCurve,
    name: str,
    operations: list[dict[str, object]],
) -> NURBSSurface:
    """Refine ``surface`` along ``axis`` until that direction *is* ``target``'s basis.

    Degree first (a patch cannot carry a seam of higher degree than its own), then
    knots: every interior breakpoint ``target`` has is inserted up to ``target``'s
    multiplicity. Both operations are exact -- the patch does not move -- and they are
    the only two that change a knot vector without changing the geometry.
    """
    current = surface
    elevation = target.degree - current.degree[axis]
    if elevation > 0:
        current = (
            current.elevate_degree(elevation, 0)
            if axis == 0
            else current.elevate_degree(0, elevation)
        )
        operations.append(
            {"patch": name, "axis": axis, "op": "elevate_degree", "times": elevation}
        )
    degree = target.degree
    count = len(current.knots[axis]) - degree - 1
    have = _axis_multiplicities(current.knots[axis], degree, count)
    want = _axis_multiplicities(target.knots, degree, len(target.control_points))
    for value in sorted(want):
        deficit = want[value] - have.get(value, 0)
        if deficit > 0:
            current = _insert_surface_knot(current, axis, value, deficit)
            operations.append(
                {
                    "patch": name,
                    "axis": axis,
                    "op": "insert_knot",
                    "value": value,
                    "times": deficit,
                }
            )
    return current


def _max_patch_deviation(first: NURBSSurface, second: NURBSSurface, samples: int = 17) -> float:
    """Largest distance between two patches at the same ``(u, v)`` on a sample grid."""
    (u_low, u_high), (v_low, v_high) = first.domain
    worst = 0.0
    for row in range(samples):
        u = u_low + (u_high - u_low) * row / (samples - 1)
        for column in range(samples):
            v = v_low + (v_high - v_low) * column / (samples - 1)
            worst = max(worst, _max_distance(first.evaluate(u, v), second.evaluate(u, v)))
    return worst


def align_patch_seam(
    first: NURBSSurface,
    first_side: str,
    second: NURBSSurface,
    second_side: str,
    *,
    reverse_second: bool = False,
    samples: int = 129,
    tolerance: float = 1e-12,
) -> SeamAlignment:
    """Make the edge ``first_side`` of ``first`` and ``second_side`` of ``second`` one
    and the same *patch data*.

    :func:`align_edges` aligns two loose curves; a stitcher cannot use those curves,
    because the patches still carry their old boundary rows and nothing in the patch
    knows about the aligned edge. This function closes that loop: it aligns the two
    boundary curves (degree, then knot vector) and then refines each patch along its
    own edge direction until the patch's knot vector equals the aligned edge's --
    by degree elevation when the seam is of higher degree, by knot insertion
    otherwise. Both are exact, so the patch geometry is untouched and only the
    representation along that direction grows.

    The axes come from the sides, not from a guess: ``"v0"``/``"v1"`` are edges in
    ``u`` (axis 0), ``"u0"``/``"u1"`` are edges in ``v`` (axis 1). ``reverse_second``
    says that the neighbour traverses the shared edge the other way, exactly as in
    :func:`boundary_curve`.

    The returned :class:`SeamAlignment` carries the two refined patches, the aligned
    edges and a report that measures the result instead of asserting it:
    ``watertight`` (equal degree, equal knots, equal control points and weights within
    ``tolerance``), ``max_gap`` at the same parameters, ``max_patch_deviation`` (how
    far each refined patch moved from its input -- round-off, if the refinement is
    exact) and ``representation_matches_target`` (whether the refinement reached the
    aligned edge's basis on both sides). Two edges that are *different* curves are
    still aligned in representation and still reported with their gap: refinement
    cannot invent geometry that is not there.
    """
    for side in (first_side, second_side):
        if side not in SIDES:
            raise ValueError(f"side must be one of {SIDES}, got {side!r}")
    first_axis = _EDGE_AXIS[first_side]
    second_axis = _EDGE_AXIS[second_side]

    first_edge = boundary_curve(first, first_side)
    second_edge = boundary_curve(second, second_side, reverse=reverse_second)
    target_first, target_second = align_edges(first_edge, second_edge)

    operations: list[dict[str, object]] = []
    refined_first = _refine_axis(first, first_axis, target_first, "first", operations)
    refined_second = _refine_axis(second, second_axis, target_second, "second", operations)

    new_first_edge = boundary_curve(refined_first, first_side)
    new_second_edge = boundary_curve(refined_second, second_side, reverse=reverse_second)

    report = edge_report(new_first_edge, new_second_edge, samples=samples, tolerance=tolerance)
    report["sides"] = [first_side, second_side]
    report["axes"] = [first_axis, second_axis]
    report["reversed_second"] = reverse_second
    report["representation_matches_target"] = [
        new_first_edge.degree == target_first.degree and new_first_edge.knots == target_first.knots,
        new_second_edge.degree == target_second.degree
        and new_second_edge.knots == target_second.knots,
    ]
    report["max_patch_deviation"] = [
        _max_patch_deviation(first, refined_first),
        _max_patch_deviation(second, refined_second),
    ]
    report["operations"] = len(operations)
    return SeamAlignment(
        first=refined_first,
        second=refined_second,
        first_edge=new_first_edge,
        second_edge=new_second_edge,
        report=report,
        operations=tuple(operations),
    )
