"""Sketch constraint solving: GAP-CAD-006 acceptance, and the AD kernel underneath it.

The acceptance criterion from ``docs/oss_inventory_cad/GAP-ANCHORS.md`` is a sentence a
machine can decide: ten fixtures including under- and over-constrained ones, convergence,
the correct degree of freedom count when under-constrained, and positions within 1e-10 of
the analytic value. Each of those clauses is a test below. The invariants are asserted, not
just executed: the null space is checked to actually annihilate the Jacobian, the reported
DOF is checked against the null-space dimension, and the residuals are recomputed at the
returned coordinates instead of trusting the solver's own report.
"""
from __future__ import annotations

import math
import random

import pytest

from tensorcad.autodiff import Dual, atan2, gradient, sqrt
from tensorcad.sketch import (
    Angle,
    Distance,
    Fix,
    Horizontal,
    InconsistentConstraintsError,
    Sketch,
    SketchError,
    Vertical,
)
from tensorcad.sketch_cases import CASES, run_all, run_case

# --- the autodiff kernel ---------------------------------------------------------------
ANALYTIC_CASES = [
    # (function, point, value, gradient)
    (lambda p: p[0] * p[1] + p[1] * p[1], (2.0, 3.0), 15.0, [3.0, 2.0 + 6.0]),
    (lambda p: sqrt(p[0] * p[0] + p[1] * p[1]), (3.0, 4.0), 5.0, [0.6, 0.8]),
    (lambda p: atan2(p[1], p[0]), (1.0, 1.0), math.pi / 4, [-0.5, 0.5]),
    (lambda p: p[0] ** 3 - 4.0 / p[1], (2.0, 4.0), 7.0, [12.0, 0.25]),
    (lambda p: 5.0 - p[0], (1.5, 0.0), 3.5, [-1.0, 0.0]),
]


@pytest.mark.parametrize("function,point,value,grad", ANALYTIC_CASES)
def test_forward_ad_matches_closed_form(function, point, value, grad) -> None:
    measured_value, measured_grad = gradient(function, point)
    assert measured_value == pytest.approx(value, rel=0, abs=1e-15)
    for got, want in zip(measured_grad, grad):
        assert got == pytest.approx(want, rel=0, abs=1e-15)


def test_constant_has_empty_gradient_and_scales_do_nothing() -> None:
    value, grad = gradient(lambda p: Dual.constant(3.0) * 2.0 + 7.0, (1.0, 2.0))
    assert value == 13.0
    assert grad == [0.0, 0.0]


def test__pow_and_reciprocal_derivatives() -> None:
    _, grad = gradient(lambda p: 1.0 / p[0], (2.0,))
    assert grad[0] == pytest.approx(-0.25, rel=1e-15)
    _, grad = gradient(lambda p: p[0] ** 3, (2.0,))
    assert grad[0] == pytest.approx(12.0, rel=1e-15)


# --- fixtures --------------------------------------------------------------------------
def test_the_fixture_set_covers_the_acceptance_clauses() -> None:
    assert len(CASES) >= 10
    statuses = {case.status for case in CASES}
    assert {"fully_constrained", "under_constrained", "inconsistent"} <= statuses
    assert any(case.redundant_rows > 0 for case in CASES)


@pytest.mark.parametrize("case", CASES, ids=[case.name for case in CASES])
def test_case_meets_its_analytic_answer(case) -> None:
    outcome = run_case(case)
    assert outcome.status == case.status
    assert outcome.dof == case.dof, f"{case.name}: reported DOF does not match the fixture"
    assert outcome.null_space_dimension == outcome.dof, "null-space basis and DOF disagree"
    assert outcome.redundant_rows >= case.redundant_rows
    if case.status == "inconsistent":
        assert outcome.raised == "InconsistentConstraintsError"
        assert outcome.max_residual > 1e-9
        assert outcome.blocking, "an inconsistent sketch must name the blocking constraints"
    else:
        assert outcome.max_residual <= 1e-9
        assert outcome.raised == ""
    # The acceptance bound: 1e-10 against the analytic positions.
    assert outcome.max_position_error <= 1e-10, f"{case.name}: {outcome.max_position_error}"
    assert outcome.max_radius_error <= 1e-10, f"{case.name}: radius {outcome.max_radius_error}"


def test_rectangles_redundancy_is_not_reported_as_failure() -> None:
    """The same dimension twice is a modelling convenience, not an over-constrained error."""
    outcome = next(o for o in run_all() if o.case == "redundant")
    assert outcome.status == "fully_constrained"
    assert outcome.dof == 0
    assert outcome.redundant_rows == 2
    assert outcome.max_residual <= 1e-12


def test_inconsistent_sketch_raises_with_the_blocking_rows() -> None:
    case = next(c for c in CASES if c.status == "inconsistent")
    sketch = case.build()
    with pytest.raises(InconsistentConstraintsError) as caught:
        sketch.solve()
    error = caught.value
    assert error.max_residual > 1e-9
    assert "pin B" in error.blocking and "|AB| = 5" in error.blocking
    assert error.result is not None and error.result.status == "inconsistent"


def test_non_strict_solve_returns_the_conflict_instead_of_raising() -> None:
    case = next(c for c in CASES if c.status == "inconsistent")
    sketch = case.build()
    result = sketch.solve(strict=False)
    assert result.status == "inconsistent"
    assert result.blocking


# --- invariants that do not depend on the fixture answers ------------------------------
@pytest.mark.parametrize("case", [c for c in CASES if c.status != "inconsistent"],
                         ids=[c.name for c in CASES if c.status != "inconsistent"])
def test_reported_null_space_really_annihilates_the_jacobian(case) -> None:
    """A degree of freedom is only reported if there is a direction the Jacobian sees as free."""
    sketch = case.build()
    result = sketch.solve()
    jacobian = sketch.jacobian()
    for vector in result.null_space:
        for row in jacobian:
            projection = sum(a * b for a, b in zip(row, vector))
            assert abs(projection) <= 1e-8, "null-space vector is not in the Jacobian null space"


def test_residuals_recomputed_at_the_returned_coordinates_are_zero() -> None:
    """Independent re-check: the solver's own residual report is not trusted here."""
    for case in CASES:
        if case.status == "inconsistent":
            continue
        sketch = case.build()
        result = sketch.solve()
        assert max(abs(v) for v in sketch.residuals(result.values)) <= 1e-9


def test_a_free_sketch_has_two_degrees_of_freedom_per_point() -> None:
    sketch = Sketch()
    sketch.point(1.0, 2.0, "A")
    sketch.point(3.0, 4.0, "B")
    result = sketch.solve()
    assert result.status == "under_constrained"
    assert result.dof == 4
    assert len(result.null_space) == 4


def test_empty_sketch_is_fully_constrained() -> None:
    result = Sketch().solve()
    assert result.status == "fully_constrained" and result.dof == 0


def test_perturbed_starts_converge() -> None:
    """25 deterministic restarts per solvable fixture: the solve is not tied to the start."""
    rng = random.Random(20261112)
    solvable = [c for c in CASES if c.status != "inconsistent"]
    attempts = 0
    for case in solvable:
        for _ in range(25):
            sketch = case.build()
            sketch.perturb(0.25, rng)
            result = sketch.solve(strict=False)
            attempts += 1
            assert result.max_residual <= 1e-9, f"{case.name}: stalled at {result.max_residual}"
            for name, expected in case.positions.items():
                solved = result.positions[name]
                assert math.hypot(solved[0] - expected[0], solved[1] - expected[1]) <= 1e-7
    assert attempts == 25 * len(solvable)


def test_rank_and_position_answers_do_not_depend_on_the_drawing_units() -> None:
    """Row-normalised rank: the same rectangle in millimetres and kilometres."""

    def rectangle(scale: float) -> tuple[Sketch, dict[str, tuple[float, float]]]:
        sketch = Sketch()
        a = sketch.point(0.0, 0.0, "A")
        b = sketch.point(4.0 * scale, 0.3 * scale, "B")
        c = sketch.point(4.2 * scale, 3.0 * scale, "C")
        d = sketch.point(0.1 * scale, 3.0 * scale, "D")
        ab, bc, cd, da = (sketch.line(a, b), sketch.line(b, c), sketch.line(c, d),
                          sketch.line(d, a))
        sketch.add(Fix(a, 0.0, 0.0), Horizontal(ab), Vertical(bc), Horizontal(cd), Vertical(da),
                   Distance(a, b, 4.0 * scale), Distance(b, c, 3.0 * scale))
        return sketch, {"A": (0.0, 0.0), "B": (4.0 * scale, 0.0),
                        "C": (4.0 * scale, 3.0 * scale), "D": (0.0, 3.0 * scale)}

    for scale in (1e-3, 1.0, 1e3):
        sketch, expected = rectangle(scale)
        # The convergence tolerance is stated in the drawing's own units, exactly as the
        # dimensions are: asking for 1e-12 m on a sketch drawn in kilometres is a
        # different request from asking for 1e-12 m on one drawn in millimetres.
        result = sketch.solve(tol=1e-12 * scale)
        assert result.dof == 0 and result.status == "fully_constrained"
        for name, want in expected.items():
            got = result.positions[name]
            error = math.hypot(got[0] - want[0], got[1] - want[1])
            assert error / (4.0 * scale) <= 1e-10, f"scale {scale}: {error}"


def test_solving_twice_is_idempotent() -> None:
    sketch = Sketch()
    a = sketch.point(0.0, 0.0, "A")
    b = sketch.point(1.0, 1.0, "B")
    ab = sketch.line(a, b, "AB")
    sketch.add(Fix(a, 0.0, 0.0), Horizontal(ab), Distance(a, b, 2.0))
    first = sketch.solve()
    second = sketch.solve()
    assert first.positions == second.positions
    assert second.max_residual <= 1e-12 and second.dof == 0


def test_initial_vector_of_the_wrong_length_is_rejected() -> None:
    sketch = Sketch()
    sketch.point(0.0, 0.0, "A")
    with pytest.raises(SketchError):
        sketch.solve(initial=[0.0, 0.0, 0.0])


def test_angle_residual_uses_degrees_and_is_signed() -> None:
    """The reported angle is the signed angle from the first line to the second."""
    sketch = Sketch()
    o = sketch.point(0.0, 0.0, "O")
    r = sketch.point(1.0, 0.0, "R")
    tip = sketch.point(1.0, -1.0, "TIP")
    reference, arm = sketch.line(o, r, "OR"), sketch.line(o, tip, "ARM")
    sketch.add(Fix(o, 0.0, 0.0), Fix(r, 1.0, 0.0), Distance(o, tip, 2.0),
               Angle(reference, arm, -45.0))
    result = sketch.solve()
    assert result.dof == 0
    assert result.positions["TIP"] == pytest.approx((math.sqrt(2.0), -math.sqrt(2.0)), abs=1e-12)