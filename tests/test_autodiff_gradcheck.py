"""Gradcheck for the forward-mode AD kernel that ``tensorcad.sketch`` builds its Jacobian on.

Why this file exists: ``docs/verification_status_m16807.md`` lists every module with an
analytic or literature anchor, and ``tensorcad.autodiff`` is the one module in the package
with *no* dedicated verification document -- it is only exercised indirectly through the
sketch fixtures. But ``docs/sketch.md`` accepts a solver whose positions land within
``1e-10`` of the analytic value, and that bound is unreachable if the derivative kernel is
wrong, so the kernel deserves its own anchor.

The oracle here is the **complex-step derivative** ``Im f(x + i h) / h`` with ``h = 1e-30``:
for a real-analytic expression it carries no subtractive cancellation, so it agrees with the
exact derivative to rounding error. That makes it an oracle independent of the AD code path
(it runs on Python ``complex``, never on :class:`Dual`), which is what a gradcheck needs --
comparing AD against a central difference would only prove both agree, while a central
difference has its own O(h^2) truncation error that can mask a wrong sign or factor at
coarse ``h``.

Every assertion below is an invariant, not an execution: values are compared to the oracle,
the behaviour at documented singularities is asserted to raise, and the quantitative claim
in the module docstring (an exact Jacobian is what makes the 1e-10 bound reachable) is
turned into a measurement of the residual floor a Newton iteration actually reaches.
"""
from __future__ import annotations

import cmath
import math

import pytest

from tensorcad.autodiff import Dual, atan2, cos, gradient, sin, sqrt, variable

# Complex-step is exact to rounding for analytic expressions; h must stay far above the
# underflow floor and far below any feature length so that f(x+ih) is still resolved.
H = 1e-30


def complex_step(function, point, index):
    """d f / d x_index from the imaginary part of a complex perturbation."""
    shift = complex(0.0, H)
    args = [complex(v) for v in point]
    args[index] += shift
    return function(args).imag / H


# --- the gradcheck table ---------------------------------------------------------------
# Each case carries the AD expression and one algebraically identical expression on Python
# complex numbers. The two are written separately on purpose: the oracle must not run
# through the code it is checking.
def _poly_dual(p):
    return p[0] ** 3 - 2.0 * p[0] * p[1] + p[1] ** 2 + 5.0


def _poly_cx(z):
    return z[0] ** 3 - 2.0 * z[0] * z[1] + z[1] ** 2 + 5.0


def _rational_dual(p):
    return (p[0] * p[0] + 1.0) / (p[1] - 0.3)


def _rational_cx(z):
    return (z[0] * z[0] + 1.0) / (z[1] - 0.3)


def _trig_dual(p):
    return sin(p[0]) * cos(p[1]) + atan2(p[1], p[0])


def _trig_cx(z):
    return cmath.sin(z[0]) * cmath.cos(z[1]) + cmath.atan(z[1] / z[0])


def _length_dual(p):
    return sqrt(p[0] * p[0] + p[1] * p[1])


def _length_cx(z):
    return cmath.sqrt(z[0] * z[0] + z[1] * z[1])


def _fractional_dual(p):
    return p[0] ** 1.5 - p[1] ** -2.0


def _fractional_cx(z):
    return z[0] ** 1.5 - z[1] ** -2.0


def _chain_dual(p):
    # sin(atan2(x, y)) * hypot(x, y) == x identically, so both the value and the gradient
    # have an answer that no finite-difference oracle could supply: (1, 0) exactly.
    return sin(atan2(p[0], p[1])) * sqrt(p[0] * p[0] + p[1] * p[1])


def _chain_cx(z):
    return cmath.sin(cmath.atan(z[0] / z[1])) * cmath.sqrt(z[0] * z[0] + z[1] * z[1])


CASES = [
    ("poly", _poly_dual, _poly_cx, (0.7, -1.3)),
    ("rational", _rational_dual, _rational_cx, (1.1, 0.95)),
    ("trig_atan2", _trig_dual, _trig_cx, (0.6, 1.45)),
    ("hypot_sqrt", _length_dual, _length_cx, (-1.7, 0.35)),
    ("fractional_pow", _fractional_dual, _fractional_cx, (2.4, -0.8)),
]


@pytest.mark.parametrize("name,dual_fn,cx_fn,point", CASES, ids=[c[0] for c in CASES])
def test_gradient_matches_complex_step(name, dual_fn, cx_fn, point):
    value, grad = gradient(dual_fn, point)

    oracle = [complex_step(cx_fn, point, i) for i in range(len(point))]
    expected_value = cx_fn([complex(v) for v in point]).real

    assert value == pytest.approx(expected_value, rel=0.0, abs=1e-13), name
    for i, (got, ref) in enumerate(zip(grad, oracle)):
        assert got == pytest.approx(ref, rel=1e-12, abs=1e-12), f"{name} d/dx{i}"
        # The bound that makes this a *gradcheck* and not a smoke test: the oracle is
        # exact to rounding, so anything looser than 1e-12 relative would hide a real bug.
        assert abs(got - ref) <= 1e-12 * max(1.0, abs(ref)), f"{name} d/dx{i} loose"


def test_chain_identity_recovers_parameter_and_unit_gradient():
    point = (0.37, -1.21)
    value, grad = gradient(_chain_dual, point)

    # An algebraic identity, so the reference is exact rather than another numerics path.
    assert value == pytest.approx(point[0], rel=0.0, abs=1e-15)
    assert grad[0] == pytest.approx(1.0, rel=0.0, abs=1e-14)
    assert grad[1] == pytest.approx(0.0, rel=0.0, abs=1e-15)


# --- structural semantics --------------------------------------------------------------
def test_independent_variables_do_not_share_gradient_entries():
    x = variable(0, 2.0)
    y = variable(1, 3.0)
    product = x * y
    assert product.value == pytest.approx(6.0)
    assert product.grad == pytest.approx({0: 3.0, 1: 2.0})


def test_constant_coercion_on_both_sides_keeps_gradients():
    x = variable(0, 4.0)
    assert (x + 1.0).grad == {0: 1.0}
    assert (1.0 + x).grad == {0: 1.0}
    assert (3.0 - x).grad == {0: -1.0}
    assert (x - 3.0).grad == {0: 1.0}
    assert (2.0 * x).grad == {0: 2.0}
    assert (8.0 / x).grad == pytest.approx({0: -0.5})
    assert (-x).grad == {0: -1.0}


def test_constant_carries_no_gradient():
    constant = Dual.constant(2.5)
    assert constant.value == 2.5
    assert constant.grad == {}
    # Adding a bare constant to an expression must not introduce spurious entries.
    assert (variable(0, 1.0) + constant).grad == {0: 1.0}


def test_gradient_is_copied_on_construction():
    source = {0: 2.0}
    dual = Dual(1.0, source)
    source[0] = 99.0
    assert dual.grad == {0: 2.0}, "caller mutation must not reach the dual's gradient"


def test_repeated_index_accumulates():
    x = variable(0, 1.5)
    squared_twice = x * x + x * x
    assert squared_twice.grad == pytest.approx({0: 4.0 * 1.5})


def test_gradient_returns_dense_vector_of_parameter_length():
    value, grad = gradient(lambda p: p[0] * p[0], [3.0, 7.0])
    assert value == pytest.approx(9.0)
    assert grad == pytest.approx([6.0, 0.0]), "unused parameter must report derivative 0"


# --- documented failure modes ----------------------------------------------------------
def test_sqrt_of_negative_raises_value_error():
    with pytest.raises(ValueError, match="negative"):
        sqrt(Dual(-1e-9))


def test_sqrt_at_zero_raises_instead_of_returning_infinity():
    with pytest.raises(ZeroDivisionError, match="sqrt"):
        sqrt(Dual(0.0))


def test_division_by_zero_raises():
    with pytest.raises(ZeroDivisionError, match="division by zero"):
        Dual(1.0) / Dual(0.0)
    with pytest.raises(ZeroDivisionError, match="division by zero"):
        1.0 / Dual(0.0)


def test_atan2_at_origin_raises():
    with pytest.raises(ZeroDivisionError, match="atan2"):
        atan2(Dual(0.0), Dual(0.0))


def test_gradient_rejects_a_plain_float_result():
    with pytest.raises(TypeError, match="must return a Dual"):
        gradient(lambda p: 1.0, [1.0])


# --- the quantitative claim: what the exact derivative actually buys --------------------
# Measured, not asserted from the docstring. Two statements are separable:
#   (a) the AD derivative is exact to rounding, a central difference is not and cannot be
#       made uniformly accurate by tuning h; (b) on the benign systems tested below, a
#       tuned central difference still drives Newton to machine zero -- so the strong form
#       of the module docstring ("a differenced Jacobian stops at the truncation error")
#       is *not* reproduced here and must not be quoted as if it were. (b) is asserted as
#       a negative result so that a future reader cannot overclaim it.
def _steep(x):
    """A residual with the shape the sketch solver sees: a length root, not a polynomial."""
    return math.sqrt((x - 1.3) * (x - 1.3) + 0.49) - 0.9


def _steep_dual(p):
    return sqrt((p[0] - 1.3) * (p[0] - 1.3) + 0.49) - 0.9


def _steep_cx(z):
    return cmath.sqrt((z[0] - 1.3) * (z[0] - 1.3) + 0.49) - 0.9


def _central_difference(function, x, h):
    return (function(x + h) - function(x - h)) / (2.0 * h)


def test_central_difference_error_floor_vs_ad():
    """AD error is 0; a central difference has a U-shaped error curve and no safe h."""
    x = 0.7
    exact = complex_step(_steep_cx, (x,), 0)
    _, grad = gradient(_steep_dual, [x])
    ad_error = abs(grad[0] - exact)

    steps = (1e-2, 1e-4, 1e-5, 1e-6, 1e-8, 1e-10)
    fd_errors = {h: abs(_central_difference(_steep, x, h) - exact) for h in steps}

    # (a) the kernel itself: exact at the rounding level of the arithmetic.
    assert ad_error == 0.0, f"AD disagreed with the complex-step oracle by {ad_error:.3e}"
    # truncation-dominated regime: error falls like h^2.
    assert fd_errors[1e-2] > fd_errors[1e-4] > fd_errors[1e-5]
    # cancellation-dominated regime: smaller h makes it *worse* again, so no h is safe.
    assert fd_errors[1e-8] > fd_errors[1e-5]
    # and the tuned best case still sits two orders above the AD error floor.
    best = min(fd_errors.values())
    assert 1e-12 < best < 1e-9, f"best central-difference error {best:.3e}"


def test_newton_with_ad_jacobian_reaches_machine_zero_quadratically():
    """Quadratic convergence and a machine-zero residual floor, using the AD derivative."""
    x = 4.0
    residuals = []
    for _ in range(12):
        value, grad = gradient(_steep_dual, [x])
        residuals.append(abs(value))
        x -= value / grad[0]

    assert residuals[-1] == 0.0 or residuals[-1] <= 1e-16, (
        f"AD-Newton stalled at {residuals[-1]:.3e}"
    )
    tail = [r for r in residuals if 0.0 < r < 1e-2]
    assert len(tail) >= 2
    ratios = [math.log10(tail[i + 1]) / math.log10(tail[i]) for i in range(len(tail) - 1)]
    assert min(ratios) >= 1.8, f"not quadratic: {ratios}"


# --- a geometry-typical 2x2 system: the Jacobian the sketch solver actually inverts -----
THETA = math.pi / 6.0
TARGET = (math.cos(THETA), math.sin(THETA))


def _r_length(p):
    return sqrt(p[0] * p[0] + p[1] * p[1]) - 1.0


def _r_angle(p):
    return atan2(p[1], p[0]) - THETA


def _r_length_f(z):
    return math.hypot(z[0], z[1]) - 1.0


def _r_angle_f(z):
    return math.atan2(z[1], z[0]) - THETA


def _r_length_cx(z):
    return cmath.sqrt(z[0] * z[0] + z[1] * z[1]) - 1.0


def _r_angle_cx(z):
    # atan2 and atan differ by a branch constant, so their derivatives agree locally --
    # which is all a derivative oracle needs (the sample point is off the branch cut).
    return cmath.atan(z[1] / z[0]) - THETA


def _jacobian_ad(point):
    return [gradient(_r_length, point)[1], gradient(_r_angle, point)[1]]


def _jacobian_fd(point, h=1e-6):
    rows = []
    for f in (_r_length_f, _r_angle_f):
        rows.append(
            [
                (f([point[0] + h, point[1]]) - f([point[0] - h, point[1]])) / (2.0 * h),
                (f([point[0], point[1] + h]) - f([point[0], point[1] - h])) / (2.0 * h),
            ]
        )
    return rows


def _solve_newton(jacobian, start=(1.6, -0.5), iterations=30):
    point = list(start)
    for _ in range(iterations):
        jac = jacobian(point)
        residual = [_r_length_f(point), _r_angle_f(point)]
        det = jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0]
        if det == 0.0:
            break
        step = [
            (-residual[0] * jac[1][1] + residual[1] * jac[0][1]) / det,
            (-jac[0][0] * residual[1] + jac[1][0] * residual[0]) / det,
        ]
        point = [point[0] + step[0], point[1] + step[1]]
    error = math.hypot(point[0] - TARGET[0], point[1] - TARGET[1])
    return point, error


def test_geometry_system_jacobian_entries_are_exact():
    """AD Jacobian vs the complex-step oracle on the length+angle residual pair."""
    point = (1.6, -0.5)
    ad = _jacobian_ad(point)
    for row, functions in zip(ad, (_r_length_cx, _r_angle_cx)):
        for column in range(2):
            oracle = complex_step(lambda z, f=functions: f(z), point, column)
            assert row[column] == pytest.approx(oracle, rel=1e-13, abs=1e-15)


def test_geometry_system_solves_to_analytic_coordinates_with_ad_jacobian():
    point, error = _solve_newton(_jacobian_ad)
    assert error <= 1e-14, f"position error {error:.3e} above the sketch acceptance bound"
    assert point[0] == pytest.approx(TARGET[0], abs=1e-13)
    assert point[1] == pytest.approx(TARGET[1], abs=1e-13)


def test_differenced_jacobian_error_is_the_same_size_as_the_acceptance_bound():
    """The honest quantitative case for AD: a *margin*, not a demonstrable FD failure.

    ``docs/sketch.md`` accepts positions within 1e-10 of the analytic value. The tuned
    central difference puts ~1e-10 of error into the Jacobian entries -- the same order as
    the bound it is supposed to police -- while the AD entries carry none. That is the
    measurement; the stronger claim that FD-Newton *cannot* converge is tested below and
    found false on these systems.
    """
    point = (1.6, -0.5)
    ad = _jacobian_ad(point)
    fd = _jacobian_fd(point)
    worst = max(abs(ad[i][j] - fd[i][j]) for i in range(2) for j in range(2))

    assert 1e-11 < worst < 1e-9, f"central-difference Jacobian error {worst:.3e}"


def test_honest_negative_result_fd_newton_also_converges_here():
    """Negative result, asserted: on these benign systems the FD Jacobian is good enough.

    Reproducing this keeps the ledger honest -- the AD kernel's justification measured
    here is derivative exactness (and therefore a margin at the 1e-10 acceptance bound),
    not a failure of a tuned central difference on well-conditioned residuals.
    """
    _, error = _solve_newton(_jacobian_fd)
    assert error <= 1e-12, (
        "reported negative result no longer holds: FD-Newton now misses the analytic "
        f"coordinates by {error:.3e} -- the AD justification needs to be restated"
    )
