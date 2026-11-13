"""Forward-mode automatic differentiation: one scalar that carries its own gradient.

The constraint solver needs exact Jacobian rows, not finite differences: a Newton step
built on a differenced Jacobian stops at the difference truncation error, and the sketch
acceptance bound (positions within ``1e-10`` of the analytic value) is below what a
``1e-6`` central difference can promise after the solver amplifies it. Differentiating
each residual by hand is possible but hand-written partials are exactly the kind of code
that silently disagrees with the residual it is supposed to differentiate.

So the residuals are written once, as ordinary arithmetic, and the derivative comes from
the arithmetic: every value is a :class:`Dual` -- a float plus a sparse gradient
``{parameter index: d value / d parameter}`` -- and ``+ - * / sqrt atan2 sin cos`` carry
the chain rule. Because the gradient is a dict keyed by parameter index, independence is
implicit: two variables are unrelated exactly when their dicts share no key, and a
constant is simply a dual with an empty gradient. No tapes, no graph, no vectorised
sweep of one tangent direction at a time -- this is forward mode over the whole gradient
vector at once, which is what a residual with a handful of variables wants.

Scope: real-valued scalar arithmetic and the handful of elementary functions the sketch
constraints use. There is no matrix layer and no second derivative; ``gradient`` is the
only entry point a caller normally needs.
"""
from __future__ import annotations

import math
from collections.abc import Callable, Mapping, Sequence

__all__ = ["Dual", "gradient", "variable"]


class Dual:
    """A real value together with its gradient with respect to named parameters.

    ``grad`` maps a parameter index to the partial derivative. It is copied on
    construction so that a caller cannot mutate a gradient through the object it passed
    in, and it stays sparse: adding a constant touches nothing.
    """

    __slots__ = ("value", "grad")

    def __init__(self, value: float, grad: Mapping[int, float] | None = None) -> None:
        self.value = float(value)
        self.grad: dict[int, float] = dict(grad) if grad else {}

    # -- construction -----------------------------------------------------------------
    @classmethod
    def constant(cls, value: float) -> "Dual":
        return cls(value)

    @classmethod
    def variable(cls, index: int, value: float) -> "Dual":
        return cls(value, {index: 1.0})

    # -- helpers ----------------------------------------------------------------------
    @staticmethod
    def _coerce(value: "Dual | float") -> "Dual":
        """A float is a constant dual; this keeps every operator one code path."""
        return value if isinstance(value, Dual) else Dual(value)

    def __repr__(self) -> str:  # pragma: no cover - debugging aid
        return f"Dual({self.value!r}, {self.grad!r})"

    def _scaled(self, factor: float) -> "Dual":
        return Dual(self.value, {k: factor * v for k, v in self.grad.items()})

    def _plus(self, other: "Dual", sign: float) -> "Dual":
        grad = dict(self.grad)
        for k, v in other.grad.items():
            grad[k] = grad.get(k, 0.0) + sign * v
        return Dual(self.value + sign * other.value, grad)

    # -- arithmetic -------------------------------------------------------------------
    def __add__(self, other: "Dual | float") -> "Dual":
        return self._plus(self._coerce(other), 1.0)

    __radd__ = __add__

    def __sub__(self, other: "Dual | float") -> "Dual":
        return self._plus(self._coerce(other), -1.0)

    def __rsub__(self, other: "Dual | float") -> "Dual":
        return self._coerce(other)._plus(self, -1.0)

    def __neg__(self) -> "Dual":
        return self._scaled(-1.0)

    def __mul__(self, other: "Dual | float") -> "Dual":
        other = self._coerce(other)
        grad: dict[int, float] = {k: v * other.value for k, v in self.grad.items()}
        for k, v in other.grad.items():
            grad[k] = grad.get(k, 0.0) + v * self.value
        return Dual(self.value * other.value, grad)

    __rmul__ = __mul__

    def __truediv__(self, other: "Dual | float") -> "Dual":
        other = self._coerce(other)
        if other.value == 0.0:
            raise ZeroDivisionError("division by zero in a differentiated expression")
        inv = 1.0 / other.value
        grad: dict[int, float] = {k: v * inv for k, v in self.grad.items()}
        scale = -self.value * inv * inv
        for k, v in other.grad.items():
            grad[k] = grad.get(k, 0.0) + scale * v
        return Dual(self.value * inv, grad)

    def __rtruediv__(self, other: "Dual | float") -> "Dual":
        return self._coerce(other).__truediv__(self)

    def __pow__(self, power: float) -> "Dual":
        value = self.value**power
        # d/dx x**p = p * x**(p-1); at x = 0 with p < 1 this is singular and the caller
        # gets a ZeroDivisionError rather than a quiet nan.
        return self._scaled(power * self.value ** (power - 1.0))._replace_value(value)

    def _replace_value(self, value: float) -> "Dual":
        return Dual(value, self.grad)

    def __float__(self) -> float:
        return self.value


def variable(index: int, value: float) -> Dual:
    """A dual whose gradient is one unit in ``index``."""
    return Dual.variable(index, value)


def sqrt(value: Dual) -> Dual:
    """``sqrt`` with its derivative ``1 / (2 sqrt(x))``; negative input raises."""
    if value.value < 0.0:
        raise ValueError(f"sqrt of a negative value ({value.value!r}) in a residual")
    root = math.sqrt(value.value)
    if root == 0.0:
        # The derivative is unbounded at the origin. A sketch that asks for the length
        # of a zero-length line has no well-posed Jacobian row; say so instead of
        # handing the solver an infinity that turns into a nan step.
        raise ZeroDivisionError("sqrt derivative is singular at zero")
    return value._scaled(0.5 / root)._replace_value(root)


def atan2(y: Dual, x: Dual) -> Dual:
    """Two-argument arctangent, differentiated through both arguments.

    d atan2(y, x) = (x dy - y dx) / (x**2 + y**2); the (0, 0) case is undefined and
    raises, which is the honest answer for a zero-length reference direction.
    """
    denominator = x.value * x.value + y.value * y.value
    if denominator == 0.0:
        raise ZeroDivisionError("atan2(0, 0) is undefined in a residual")
    grad: dict[int, float] = {}
    for k, v in y.grad.items():
        grad[k] = grad.get(k, 0.0) + x.value * v / denominator
    for k, v in x.grad.items():
        grad[k] = grad.get(k, 0.0) - y.value * v / denominator
    return Dual(math.atan2(y.value, x.value), grad)


def sin(value: Dual) -> Dual:
    return Dual(math.sin(value.value), {k: math.cos(value.value) * v for k, v in value.grad.items()})


def cos(value: Dual) -> Dual:
    return Dual(math.cos(value.value), {k: -math.sin(value.value) * v for k, v in value.grad.items()})


def gradient(
    function: Callable[[Sequence[Dual]], Dual], point: Sequence[float]
) -> tuple[float, list[float]]:
    """Value and dense gradient of a scalar function of a parameter vector.

    This is the entry point a caller in doubt reaches for: it turns any function written
    in ordinary arithmetic into a value and a gradient, so a residual can be checked
    against a finite difference in a test without duplicating its algebra.
    """
    params = [Dual.variable(i, v) for i, v in enumerate(point)]
    result = function(params)
    if not isinstance(result, Dual):  # a function that returns a plain float has no gradient
        raise TypeError("function must return a Dual; wrap constants with Dual(...)")
    dense = [0.0] * len(point)
    for index, value in result.grad.items():
        dense[index] = value
    return result.value, dense