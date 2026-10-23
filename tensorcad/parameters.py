"""Parameter definitions, units and validation.

Parameters are the reproducibility boundary of TensorCAD: a geometry result is
only reproducible if the numbers that produced it are known to be well typed, in
range, and stated in a unit system. Validation therefore reports *all* problems
at once (an agent fixing a model gets one actionable list rather than one error
per round trip) and refuses to guess what a bad value meant.
"""
from __future__ import annotations

import math
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from typing import Any

FLOAT = "float"
INT = "int"


@dataclass(frozen=True)
class ParameterSpec:
    """One named parameter: type, inclusive bounds, and the unit it is stated in.

    ``minimum``/``maximum`` of ``None`` mean unbounded in that direction. A
    positive *finite* bound is what stops an accidental ``1e12`` from being
    turned into an unbounded amount of memory by the generator.
    """

    name: str
    kind: str = FLOAT
    minimum: float | None = 0.0
    maximum: float | None = None
    unit: str = "m"
    description: str = ""
    exclusive_minimum: bool = False

    def check(self, value: Any) -> str | None:
        if self.kind not in (FLOAT, INT):  # guard against a typo in a spec table
            raise ValueError(f"unknown kind {self.kind!r} for parameter {self.name!r}")
        if self.kind == INT:
            # A count is never silently coerced from 11.0 to 11: a float count is
            # almost always a division that went wrong upstream, and accepting it
            # would hide that while producing a mesh nobody asked for.
            if isinstance(value, bool) or not isinstance(value, int):
                return (
                    f"{self.name} must be a plain int, got {value!r} "
                    f"({type(value).__name__}); convert explicitly with int(...)"
                )
            return self._range(value)
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            return f"{self.name} must be a real number, got {value!r}"
        if not math.isfinite(value):
            return f"{self.name} must be finite, got {value!r}"
        return self._range(float(value))

    def _range(self, value: float) -> str | None:
        if self.minimum is not None:
            if self.exclusive_minimum:
                if value <= self.minimum:
                    return f"{self.name} must be > {self.minimum:g} {self.unit}, got {value:g}"
            elif value < self.minimum:
                return f"{self.name} must be >= {self.minimum:g} {self.unit}, got {value:g}"
        if self.maximum is not None and value > self.maximum:
            return f"{self.name} must be <= {self.maximum:g} {self.unit}, got {value:g}"
        return None


def validate(specs: Sequence[ParameterSpec], values: Mapping[str, Any]) -> list[str]:
    """Return one error string per invalid value; empty means all good."""
    errors: list[str] = []
    for spec in specs:
        if spec.name not in values:
            errors.append(f"missing parameter {spec.name!r}")
            continue
        problem = spec.check(values[spec.name])
        if problem:
            errors.append(problem)
    unknown = set(values) - {spec.name for spec in specs}
    if unknown:
        errors.append(f"unknown parameters: {', '.join(sorted(unknown))}")
    return errors


def require(specs: Sequence[ParameterSpec], values: Mapping[str, Any]) -> None:
    """Raise :class:`ValueError` listing every validation problem found."""
    errors = validate(specs, values)
    if errors:
        raise ValueError("; ".join(errors))


@dataclass(frozen=True)
class RatioSpec:
    """A dimensionless ratio check: ``numerator / denominator`` must lie in range.

    These never fail a build -- a needle hull is still valid geometry -- but they
    flag numbers that no engineer meant to type, which is exactly what a research
    agent needs back as feedback.
    """

    name: str
    numerator: str
    denominator: str
    minimum: float
    maximum: float
    guidance: str = ""

    def check(self, values: Mapping[str, float]) -> str | None:
        denominator = values[self.denominator]
        if not math.isfinite(denominator) or denominator <= 0:
            return None  # the plain bounds check already reported this
        ratio = values[self.numerator] / denominator
        if self.minimum <= ratio <= self.maximum:
            return None
        hint = f" ({self.guidance})" if self.guidance else ""
        return (
            f"{self.name} = {self.numerator}/{self.denominator} = {ratio:.4g} is outside "
            f"[{self.minimum:g}, {self.maximum:g}]{hint}"
        )


def check_ratios(specs: Iterable[RatioSpec], values: Mapping[str, float]) -> list[str]:
    return [problem for problem in (spec.check(values) for spec in specs) if problem]
