"""Degenerate starts and the status contract of the sketch solver (m16386 review).

Two defects found by probing :mod:`tensorcad.sketch` against its own published
contract, and pinned here so they cannot come back:

1. **A start on the singular set of a length row was fatal.** ``Distance(q, p, 5)`` with
   ``q`` and ``p`` both at the origin -- the most natural way to write a sketch
   programmatically, and the state every sketch is in before it is dragged -- raised a
   bare ``ZeroDivisionError`` from the automatic-differentiation kernel (the derivative
   of ``|q - p|`` is ``d / |d|``, undefined at ``d = 0``). The geometry is well posed;
   only the *start* is not. The solver now steps deterministically off the singular point
   once and retries, and reports that it did (``nudged_from_singular_start``).
   A genuinely ill-posed sketch -- a direction constraint on a line whose two endpoints
   the constraints hold at the same place -- still fails, but with a :class:`SketchError`
   naming the constraint instead of an arithmetic exception.

2. **``fully_constrained`` did not imply "solved".** The status was the structural
   verdict (``dof == 0``) with no reference to whether the residual had reached ``tol``,
   so a solve that stalled between ``tol`` and ``residual_tol`` -- a band a caller widens
   by passing ``residual_tol`` -- returned ``status == "fully_constrained"`` with
   ``converged is False`` and coordinates that miss their dimensions. The verdict is now
   withheld: that case is ``"not_converged"``, its blocking rows are named, and
   ``fully_constrained`` can only be true when the residual is inside ``tol``.
"""
from __future__ import annotations

import math

import pytest

from tensorcad.sketch import (
    SINGULAR_START_NUDGE,
    Distance,
    Fix,
    Horizontal,
    Parallel,
    Sketch,
    SketchError,
)


def _distance(first: tuple[float, float], second: tuple[float, float]) -> float:
    return math.dist(first, second)


class TestASingularStartIsNotFatal:
    """The start is on the singular set; the geometry is not."""

    def test_a_distance_between_points_that_start_together_is_solved(self) -> None:
        sketch = Sketch()
        origin = sketch.point(0.0, 0.0)
        free = sketch.point(0.0, 0.0)  # exactly on top of ``origin``: |d| = 0
        sketch.add(Fix(origin, 0.0, 0.0), Distance(free, origin, 5.0))

        result = sketch.solve()

        assert result.nudged_from_singular_start is True
        assert result.converged and result.status == "under_constrained"
        assert _distance(result.positions["P0"], result.positions["P1"]) == pytest.approx(5.0, abs=1e-12)

    def test_a_triangle_from_three_coincident_points_recovers_its_side_lengths(self) -> None:
        sketch = Sketch()
        a = sketch.point(0.0, 0.0)
        b = sketch.point(0.0, 0.0)
        c = sketch.point(0.0, 0.0)
        sketch.add(Fix(a, 0.0, 0.0), Distance(b, a, 4.0), Distance(c, a, 3.0), Distance(c, b, 2.0))

        result = sketch.solve()
        positions = result.positions

        assert result.nudged_from_singular_start is True
        assert _distance(positions["P0"], positions["P1"]) == pytest.approx(4.0, abs=1e-12)
        assert _distance(positions["P0"], positions["P2"]) == pytest.approx(3.0, abs=1e-12)
        assert _distance(positions["P1"], positions["P2"]) == pytest.approx(2.0, abs=1e-12)
        # dof 1 is the triangle spinning about the pinned vertex, not an error: three
        # distances fix the shape, one vertex is pinned, the orientation about it is free.
        assert result.dof == 1
        # independent re-check: the returned vector fed back through the residuals.
        recheck = sketch.residuals(result.values)
        assert max(abs(value) for value in recheck) < 1e-12

    def test_the_escape_is_deterministic(self) -> None:
        def solve_once() -> tuple[float, ...]:
            sketch = Sketch()
            a = sketch.point(0.0, 0.0)
            b = sketch.point(0.0, 0.0)
            sketch.add(Fix(a, 0.0, 0.0), Distance(b, a, 5.0))
            return sketch.solve().values

        assert solve_once() == solve_once()

    def test_a_well_posed_start_is_not_nudged(self) -> None:
        sketch = Sketch()
        a = sketch.point(0.0, 0.0)
        b = sketch.point(1.0, 0.5)
        sketch.add(Fix(a, 0.0, 0.0), Horizontal(sketch.line(a, b)), Distance(b, a, 5.0))

        result = sketch.solve()

        assert result.nudged_from_singular_start is False
        # horizontal fixes y, the distance fixes |b - a|: the answer is locally unique
        assert result.status == "fully_constrained" and result.dof == 0
        assert result.positions["P1"] == pytest.approx((5.0, 0.0), abs=1e-12)

    def test_the_nudge_is_smaller_than_any_dimension_a_drawing_would_state(self) -> None:
        # document the magnitude the docstring claims: micrometres, not millimetres
        assert 0.0 < SINGULAR_START_NUDGE <= 1e-5

    def test_a_line_the_constraints_hold_at_zero_length_is_never_reported_as_solved(self) -> None:
        # A direction constraint on a line whose endpoints the constraints pin together is
        # ill-posed geometry. Whatever route the solver takes, it must not answer
        # "fully_constrained", and the automatic-differentiation singularity must not
        # escape as an arithmetic exception.
        sketch = Sketch()
        a = sketch.point(0.0, 0.0)
        b = sketch.point(0.0, 0.0)
        c = sketch.point(0.0, 1.0)
        d = sketch.point(1.0, 1.0)
        sketch.add(Fix(a, 0.0, 0.0), Fix(b, 0.0, 0.0), Fix(c, 0.0, 1.0), Fix(d, 1.0, 1.0))
        sketch.add(Parallel(sketch.line(a, b), sketch.line(c, d)))

        try:
            result = sketch.solve(strict=False)
        except SketchError as error:  # named constraint, not a ZeroDivisionError
            assert "parallel" in str(error)
            return
        assert result.status in {"inconsistent", "not_converged"}
        assert result.converged is False
        assert result.fully_constrained is False

    def test_the_public_residuals_report_a_singular_configuration_as_a_sketch_error(self) -> None:
        sketch = Sketch()
        a = sketch.point(0.0, 0.0)
        b = sketch.point(0.0, 0.0)
        c = sketch.point(0.0, 1.0)
        d = sketch.point(1.0, 1.0)
        sketch.add(Parallel(sketch.line(a, b), sketch.line(c, d)))

        with pytest.raises(SketchError):
            sketch.residuals()
        with pytest.raises(SketchError):
            sketch.jacobian()


class TestTheStatusIsASolutionVerdict:
    """``fully_constrained`` must never be reported for coordinates that miss a dimension."""

    @staticmethod
    def _unreachable_dimension(**solve_kwargs) -> object:
        # (3, 0) pinned and asked to be 5 m from the pinned origin: no solution exists.
        sketch = Sketch()
        origin = sketch.point(0.0, 0.0)
        point = sketch.point(3.0, 0.0)
        sketch.add(Fix(origin, 0.0, 0.0), Fix(point, 3.0, 0.0), Distance(point, origin, 5.0))
        return sketch.solve(strict=False, **solve_kwargs)

    def test_a_stalled_solve_inside_the_inconsistency_band_is_not_fully_constrained(self) -> None:
        result = self._unreachable_dimension(tol=1e-3, residual_tol=1.0)

        assert result.status == "not_converged"
        assert result.not_converged is True
        assert result.fully_constrained is False
        assert result.converged is False
        assert result.solved is False
        assert result.max_residual > 1e-3
        # the rows that are still out are named, so the caller can act on them
        assert "distance" in result.blocking

    def test_a_solve_stopped_above_tolerance_at_the_defaults_is_not_fully_constrained(self) -> None:
        # tol below the achievable round-off: the residual stops at ~1e-16 and the
        # structural verdict must not be handed out.
        sketch = Sketch()
        a = sketch.point(0.0, 0.0)
        b = sketch.point(0.0, 0.0)
        c = sketch.point(0.0, 0.0)
        sketch.add(Fix(a, 0.0, 0.0), Distance(b, a, 4.0), Distance(c, a, 3.0), Distance(c, b, 2.0))

        result = sketch.solve(tol=1e-30, residual_tol=1e3, max_iterations=1)

        assert result.converged is False
        assert result.status == "not_converged"
        assert result.fully_constrained is False

    def test_inconsistency_still_takes_precedence(self) -> None:
        result = self._unreachable_dimension()

        assert result.status == "inconsistent"
        assert "distance" in result.blocking

    def test_a_converged_solve_keeps_its_structural_verdict(self) -> None:
        sketch = Sketch()
        a = sketch.point(0.0, 0.0)
        b = sketch.point(3.0, 0.0)
        sketch.add(Fix(a, 0.0, 0.0), Fix(b, 3.0, 0.0))

        result = sketch.solve()

        assert result.status == "fully_constrained"
        assert result.converged is True
        assert result.solved is True
        assert result.blocking == ()


class TestTheTwoDefectsCompose:
    """A singular start on a sketch that cannot be satisfied: the status stays honest."""

    def test_a_stalled_solve_from_a_coincident_start_is_not_converged(self) -> None:
        sketch = Sketch()
        origin = sketch.point(0.0, 0.0)
        point = sketch.point(0.0, 0.0)
        sketch.add(Fix(origin, 0.0, 0.0), Fix(point, 0.0, 0.0), Distance(point, origin, 5.0))

        result = sketch.solve(strict=False, tol=1e-3, residual_tol=10.0)

        assert result.nudged_from_singular_start is True
        assert result.status == "not_converged"
        assert result.fully_constrained is False
