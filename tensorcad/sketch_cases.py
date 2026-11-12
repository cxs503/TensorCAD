"""Twelve sketch fixtures with analytic solutions, shared by the tests and the example.

``GAP-CAD-006``'s acceptance criterion is stated as a machine-decidable sentence: *ten
fixtures, including under-constrained and over-constrained ones, 100 % convergence, the
correct degree-of-freedom count when under-constrained, and positions within 1e-10 of the
analytic value*. This module is that fixture set, written once so the test asserts against
it and the example reports from it -- the numbers in ``docs/sketch_constraints.md`` are
produced by running these cases, never typed in.

Each case is a *fresh* sketch (builders share no state, so a test can run them in any
order) plus the analytic answer, computed from closed forms in the builder's own terms:
the SSS triangle's apex from the law of cosines, the inscribed circle between two parallel
lines from ``r = h / 2``, the 30-degree arm from ``(L cos 30, L sin 30)``. A case whose
solution is not unique (a point on a circle) states its DOF and leaves the position table
empty, because "the answer" there is the constraint residual and the remaining freedom,
not a coordinate.

Where a sketch has a discrete mirror solution the starting coordinates pick the branch;
that is a documented property of a local solver, and each case records it.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable, Mapping

from .sketch import (
    Angle,
    Distance,
    Fix,
    FixX,
    FixY,
    Horizontal,
    InconsistentConstraintsError,
    LineCircleTangency,
    Perpendicular,
    PointLineDistance,
    PointOnCircle,
    PointOnLine,
    Radius,
    Sketch,
    SolveResult,
    Vertical,
)

__all__ = ["CASES", "CaseOutcome", "SketchCase", "run_all", "run_case"]


@dataclass(frozen=True)
class SketchCase:
    """One fixture: how to build it, and what the analytic answer is."""

    name: str
    description: str
    build: Callable[[], Sketch]
    status: str  # "fully_constrained" | "under_constrained" | "inconsistent"
    dof: int
    positions: Mapping[str, tuple[float, float]] = field(default_factory=dict)
    radii: Mapping[str, float] = field(default_factory=dict)
    redundant_rows: int = 0
    note: str = ""


@dataclass(frozen=True)
class CaseOutcome:
    """What the solver did on a case, with the measured error against the analytic answer."""

    case: str
    status: str
    dof: int
    rank: int
    parameters: int
    residual_rows: int
    redundant_rows: int
    iterations: int
    max_residual: float
    max_position_error: float
    max_radius_error: float
    null_space_dimension: int
    blocking: tuple[str, ...]
    raised: str
    solved_positions: dict[str, tuple[float, float]]
    analytic_positions: dict[str, tuple[float, float]]
    descriptions: str


def _new() -> Sketch:
    return Sketch()


# -- case builders ---------------------------------------------------------------------
def _fixed_point() -> Sketch:
    """One point, pinned: the smallest sketch with a unique answer."""
    sketch = _new()
    a = sketch.point(1.0, -2.0, "A")
    sketch.add(Fix(a, 2.5, -1.25, name="lock A"))
    return sketch


def _distance_circle() -> Sketch:
    """A point held 5 m from a pinned point: one equation, two parameters left -> DOF 1."""
    sketch = _new()
    a = sketch.point(0.0, 0.0, "A")
    b = sketch.point(3.0, 0.5, "B")
    sketch.add(Fix(a, 0.0, 0.0, name="lock A"), Distance(a, b, 5.0, name="|AB| = 5"))
    return sketch


def _horizontal_arm() -> Sketch:
    """Horizontal 3 m arm off a pinned point."""
    sketch = _new()
    a = sketch.point(0.0, 0.0, "A")
    b = sketch.point(2.0, 1.0, "B")
    ab = sketch.line(a, b, "AB")
    sketch.add(
        Fix(a, 0.0, 0.0, name="lock A"),
        Horizontal(ab, name="AB horizontal"),
        Distance(a, b, 3.0, name="|AB| = 3"),
    )
    return sketch


# SSS triangle: base A-B is 4 m, the apex is 3 m from A and 2 m from B. Law of cosines
# gives the apex abscissa; the ordinate follows from the A radius. Mirror branch: +y.
_SSS_BASE = 4.0
_SSS_APEX_TO_A = 3.0
_SSS_APEX_TO_B = 2.0
_SSS_X = (_SSS_BASE**2 + _SSS_APEX_TO_A**2 - _SSS_APEX_TO_B**2) / (2.0 * _SSS_BASE)
_SSS_Y = math.sqrt(_SSS_APEX_TO_A**2 - _SSS_X**2)


def _sss_triangle() -> Sketch:
    sketch = _new()
    a = sketch.point(0.0, 0.0, "A")
    b = sketch.point(_SSS_BASE, 0.0, "B")
    c = sketch.point(2.0, 2.0, "C")
    sketch.add(
        Fix(a, 0.0, 0.0, name="lock A"),
        Fix(b, _SSS_BASE, 0.0, name="lock B"),
        Distance(a, c, _SSS_APEX_TO_A, name="|AC| = 3"),
        Distance(b, c, _SSS_APEX_TO_B, name="|BC| = 2"),
    )
    return sketch


def _rectangle() -> Sketch:
    """4 x 3 rectangle: horizontal/vertical pairs plus two dimensions, DOF 0."""
    sketch = _new()
    a = sketch.point(0.0, 0.0, "A")
    b = sketch.point(4.0, 0.3, "B")
    c = sketch.point(4.2, 3.0, "C")
    d = sketch.point(0.1, 3.0, "D")
    ab, bc, cd, da = (sketch.line(a, b, "AB"), sketch.line(b, c, "BC"),
                      sketch.line(c, d, "CD"), sketch.line(d, a, "DA"))
    sketch.add(
        Fix(a, 0.0, 0.0, name="lock A"),
        Horizontal(ab, name="AB horizontal"),
        Vertical(bc, name="BC vertical"),
        Horizontal(cd, name="CD horizontal"),
        Vertical(da, name="DA vertical"),
        Distance(a, b, 4.0, name="width = 4"),
        Distance(b, c, 3.0, name="height = 3"),
    )
    return sketch


_ARM_LENGTH = 2.0
_ARM_DEGREES = 30.0
_ARM_X = _ARM_LENGTH * math.cos(math.radians(_ARM_DEGREES))
_ARM_Y = _ARM_LENGTH * math.sin(math.radians(_ARM_DEGREES))


def _angle_arm() -> Sketch:
    """2 m arm at 30 degrees to a pinned x reference line."""
    sketch = _new()
    o = sketch.point(0.0, 0.0, "O")
    r = sketch.point(1.0, 0.0, "R")
    tip = sketch.point(1.0, 1.5, "TIP")
    reference = sketch.line(o, r, "OR")
    arm = sketch.line(o, tip, "ARM")
    sketch.add(
        Fix(o, 0.0, 0.0, name="lock O"),
        Fix(r, 1.0, 0.0, name="lock R"),
        Distance(o, tip, _ARM_LENGTH, name="arm = 2"),
        Angle(reference, arm, _ARM_DEGREES, name="arm = 30 deg"),
    )
    return sketch


def _perpendicular_chain() -> Sketch:
    """Vertical 3 m arm, then a 4 m arm perpendicular to it."""
    sketch = _new()
    a = sketch.point(0.0, 0.0, "A")
    b = sketch.point(0.4, 2.8, "B")
    c = sketch.point(5.0, 3.2, "C")
    ab, bc = sketch.line(a, b, "AB"), sketch.line(b, c, "BC")
    sketch.add(
        Fix(a, 0.0, 0.0, name="lock A"),
        Vertical(ab, name="AB vertical"),
        Distance(a, b, 3.0, name="|AB| = 3"),
        Perpendicular(ab, bc, name="AB perpendicular BC"),
        Distance(b, c, 4.0, name="|BC| = 4"),
    )
    return sketch


def _point_on_line() -> Sketch:
    """2.5 m along the 3-4-5 direction: the point on a line at a given distance."""
    sketch = _new()
    a = sketch.point(0.0, 0.0, "A")
    b = sketch.point(4.0, 3.0, "B")
    p = sketch.point(1.0, 1.0, "P")
    ab = sketch.line(a, b, "AB")
    sketch.add(
        Fix(a, 0.0, 0.0, name="lock A"),
        Fix(b, 4.0, 3.0, name="lock B"),
        PointOnLine(p, ab, name="P on AB"),
        Distance(a, p, 2.5, name="|AP| = 2.5"),
    )
    return sketch


def _point_line_offset() -> Sketch:
    """Signed 1 m offset from the 3-4-5 line, plus an abscissa anchor."""
    sketch = _new()
    a = sketch.point(0.0, 0.0, "A")
    b = sketch.point(4.0, 3.0, "B")
    p = sketch.point(1.0, 2.0, "P")
    ab = sketch.line(a, b, "AB")
    # The line direction is (4, 3) / 5, so the left unit normal is (-0.6, 0.8): a point at
    # signed distance +1 sits 1 m along that normal from the line.
    sketch.add(
        Fix(a, 0.0, 0.0, name="lock A"),
        Fix(b, 4.0, 3.0, name="lock B"),
        FixX(p, -0.6, name="P x = -0.6"),
        PointLineDistance(p, ab, 1.0, name="offset = 1"),
    )
    return sketch


def _circle_between_lines() -> Sketch:
    """Circle tangent to two parallel lines 2 m apart: radius is h / 2, centre is at h / 2."""
    sketch = _new()
    a = sketch.point(0.0, 0.0, "A")
    b = sketch.point(4.0, 0.0, "B")
    c = sketch.point(0.0, 2.0, "C")
    d = sketch.point(4.0, 2.0, "D")
    o = sketch.point(1.0, 1.3, "O")
    bottom = sketch.line(a, b, "bottom")
    top = sketch.line(c, d, "top")
    circle = sketch.circle(o, 0.7, "C")
    sketch.add(
        Fix(a, 0.0, 0.0, name="lock A"),
        Fix(b, 4.0, 0.0, name="lock B"),
        Fix(c, 0.0, 2.0, name="lock C"),
        Fix(d, 4.0, 2.0, name="lock D"),
        FixX(o, 1.0, name="centre x = 1"),
        LineCircleTangency(bottom, circle, side="left", name="tangent bottom"),
        LineCircleTangency(top, circle, side="right", name="tangent top"),
    )
    return sketch


def _point_on_circle() -> Sketch:
    """A point on a circle of radius 2 about the origin, plus a y anchor."""
    sketch = _new()
    o = sketch.point(0.0, 0.0, "O")
    p = sketch.point(1.5, 1.0, "P")
    circle = sketch.circle(o, 2.0, "C")
    sketch.add(
        Fix(o, 0.0, 0.0, name="lock O"),
        Radius(circle, 2.0, name="r = 2"),
        PointOnCircle(p, circle, name="P on circle"),
        FixY(p, 1.0, name="P y = 1"),
    )
    return sketch


def _under_constrained_line() -> Sketch:
    """A bare horizontal line: 4 parameters, 1 equation -> DOF 3."""
    sketch = _new()
    a = sketch.point(0.0, 0.0, "A")
    b = sketch.point(3.0, 0.4, "B")
    ab = sketch.line(a, b, "AB")
    sketch.add(Horizontal(ab, name="AB horizontal"))
    return sketch


def _inconsistent() -> Sketch:
    """A point pinned at (3, 0) and told to be 5 m from the origin: no solution."""
    sketch = _new()
    a = sketch.point(0.0, 0.0, "A")
    b = sketch.point(2.0, 0.0, "B")
    sketch.add(
        Fix(a, 0.0, 0.0, name="lock A"),
        Distance(a, b, 5.0, name="|AB| = 5"),
        Fix(b, 3.0, 0.0, name="pin B"),
    )
    return sketch


def _redundant() -> Sketch:
    """The same dimension stated twice: consistent, but two rows beyond the rank."""
    sketch = _new()
    a = sketch.point(0.0, 0.0, "A")
    b = sketch.point(3.2, 3.8, "B")
    sketch.add(
        Fix(a, 0.0, 0.0, name="lock A"),
        Fix(b, 3.0, 4.0, name="lock B"),
        Distance(a, b, 5.0, name="dim 1"),
        Distance(a, b, 5.0, name="dim 2"),
    )
    return sketch


CASES: tuple[SketchCase, ...] = (
    SketchCase(
        name="fixed_point",
        description="one point pinned at (2.5, -1.25)",
        build=_fixed_point,
        status="fully_constrained",
        dof=0,
        positions={"A": (2.5, -1.25)},
    ),
    SketchCase(
        name="distance_circle",
        description="point held 5 m from a pinned point (one equation, two parameters)",
        build=_distance_circle,
        status="under_constrained",
        dof=1,
        note="the answer is the 5 m residual plus one tangential motion, not a coordinate",
    ),
    SketchCase(
        name="horizontal_arm",
        description="horizontal 3 m arm off a pinned point",
        build=_horizontal_arm,
        status="fully_constrained",
        dof=0,
        positions={"A": (0.0, 0.0), "B": (3.0, 0.0)},
        note="the mirror branch B = (-3, 0) is excluded by the starting guess",
    ),
    SketchCase(
        name="sss_triangle",
        description="triangle from three side lengths (4, 3, 2 m), apex from the law of cosines",
        build=_sss_triangle,
        status="fully_constrained",
        dof=0,
        positions={"A": (0.0, 0.0), "B": (_SSS_BASE, 0.0), "C": (_SSS_X, _SSS_Y)},
        note="upper apex branch; the lower mirror is a different start",
    ),
    SketchCase(
        name="rectangle",
        description="4 x 3 rectangle from horizontal/vertical pairs and two dimensions",
        build=_rectangle,
        status="fully_constrained",
        dof=0,
        positions={"A": (0.0, 0.0), "B": (4.0, 0.0), "C": (4.0, 3.0), "D": (0.0, 3.0)},
    ),
    SketchCase(
        name="angle_arm",
        description="2 m arm at 30 degrees to a pinned x reference line",
        build=_angle_arm,
        status="fully_constrained",
        dof=0,
        positions={"O": (0.0, 0.0), "R": (1.0, 0.0), "TIP": (_ARM_X, _ARM_Y)},
    ),
    SketchCase(
        name="perpendicular_chain",
        description="vertical 3 m arm followed by a 4 m arm perpendicular to it",
        build=_perpendicular_chain,
        status="fully_constrained",
        dof=0,
        positions={"A": (0.0, 0.0), "B": (0.0, 3.0), "C": (4.0, 3.0)},
    ),
    SketchCase(
        name="point_on_line",
        description="point on the 3-4-5 line at 2.5 m from its start",
        build=_point_on_line,
        status="fully_constrained",
        dof=0,
        positions={"A": (0.0, 0.0), "B": (4.0, 3.0), "P": (2.0, 1.5)},
    ),
    SketchCase(
        name="point_line_offset",
        description="point 1 m to the left of the 3-4-5 line, x anchored",
        build=_point_line_offset,
        status="fully_constrained",
        dof=0,
        positions={"A": (0.0, 0.0), "B": (4.0, 3.0), "P": (-0.6, 0.8)},
    ),
    SketchCase(
        name="circle_between_lines",
        description="circle tangent to two parallel lines 2 m apart (r = h / 2)",
        build=_circle_between_lines,
        status="fully_constrained",
        dof=0,
        positions={"O": (1.0, 1.0)},
        radii={"C": 1.0},
    ),
    SketchCase(
        name="point_on_circle",
        description="point on a radius-2 circle with its ordinate anchored",
        build=_point_on_circle,
        status="fully_constrained",
        dof=0,
        positions={"P": (math.sqrt(3.0), 1.0)},
        radii={"C": 2.0},
    ),
    SketchCase(
        name="under_constrained_line",
        description="a bare horizontal line: 4 parameters, one equation",
        build=_under_constrained_line,
        status="under_constrained",
        dof=3,
        note="the three motions are the two translations and the endpoint split",
    ),
    SketchCase(
        name="inconsistent",
        description="point pinned at (3, 0) and asked to be 5 m from the origin",
        build=_inconsistent,
        status="inconsistent",
        dof=0,
        note="over-constrained and conflicting: the solver names both blocking rows",
    ),
    SketchCase(
        name="redundant",
        description="the same 5 m dimension stated twice over two pinned points",
        build=_redundant,
        status="fully_constrained",
        dof=0,
        positions={"A": (0.0, 0.0), "B": (3.0, 4.0)},
        redundant_rows=2,
        note="redundant but consistent: reported as redundancy, not as failure",
    ),
)


def run_case(case: SketchCase) -> CaseOutcome:
    """Build, solve, and measure one case; no comparison happens here (the caller judges)."""
    sketch = case.build()
    raised = ""
    try:
        result: SolveResult = sketch.solve()
    except InconsistentConstraintsError as error:
        raised = type(error).__name__
        assert error.result is not None
        result = error.result
    position_error = 0.0
    for name, expected in case.positions.items():
        solved = result.positions[name]
        position_error = max(position_error, math.hypot(solved[0] - expected[0],
                                                       solved[1] - expected[1]))
    radius_error = 0.0
    for name, expected in case.radii.items():
        radius_error = max(radius_error, abs(result.radii[name] - expected))
    return CaseOutcome(
        case=case.name,
        status=result.status,
        dof=result.dof,
        rank=result.rank,
        parameters=result.parameters,
        residual_rows=result.residual_rows,
        redundant_rows=result.redundant_rows,
        iterations=result.iterations,
        max_residual=result.max_residual,
        max_position_error=position_error,
        max_radius_error=radius_error,
        null_space_dimension=len(result.null_space),
        blocking=result.blocking,
        raised=raised,
        solved_positions={name: result.positions[name] for name in case.positions},
        analytic_positions=dict(case.positions),
        descriptions=case.description,
    )


def run_all() -> list[CaseOutcome]:
    return [run_case(case) for case in CASES]