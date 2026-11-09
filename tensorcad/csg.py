"""Boolean (CSG) composition of geometry fields: union, intersection, difference.

The field protocol in :mod:`tensorcad.geometry_field` answers *one* surface at a
time. Every real model is several: a hull plus a sail plus four fins, a tank
opening cut out of a compartment, the part of a bulb above the design waterline.
Until now the only way to say that in this package was to hand-write a class --
``SuboffFullField`` is exactly that, a hard-coded ``min`` of hull and appendages --
which is why "no boolean operations" has stood on the capability ledger since the
0.4.0 audit (``backlog_tensorcad_section.md``).

This module supplies the operator, once, in the standard ``min``/``max`` form:

======================  =====================================  ===========================
operation               field                                  zero set
======================  =====================================  ===========================
``union(a, b)``         ``min(a, b)``                          ``A`` union ``B``
``intersection(a, b)``  ``max(a, b)``                          ``A`` intersect ``B``
``difference(a, b)``    ``max(a, -b)``                         ``A`` minus ``B``
======================  =====================================  ===========================

**What is exact and what is not.** The *sign* of each expression is the sign of the
Boolean membership everywhere, including on the seam and in the concave creases where
two bodies meet -- ``min(a, b) < 0`` is ``a < 0 or b < 0`` for any real ``a``, ``b``.
So the zero level set, containment, mesh extraction and voxel occupancy are exact
statements about the composed solid. The *magnitude* is not a distance everywhere:
``|min(a, b)|`` is the true distance when the nearer surface belongs to the operand
that attains the minimum, and an over-estimate otherwise (a point inside ``A`` and
below ``B`` reports its distance to ``A``'s surface, not to the nearest point of the
union's surface). ``max`` is the mirror case. That is a fact about the construction,
not a defect of this implementation -- an exact CSG distance needs a nearest-point
query -- so a composed field declares ``quantity = "implicit_function"``, the same
honesty ``EllipsoidField`` and ``SuboffFullField`` already apply.

The consequence that matters for consumers is the safe direction: a ``min``/``max``
composition never *understates* the distance to the composed boundary
(``|value| >= true distance``), so marching tetrahedra, a voxel grid and a
distance-band test all stay conservative. ``tests/test_csg.py`` pins that inequality
against a brute-force nearest-surface search, and pins the closed-form volumes of
the solids the operators produce.

Boolean of a learned ``.ngo`` field is allowed and means what it says: the composed
surface is where the learned field zero-crosses. Nothing here assumes an analytic
operand -- that is the point of the protocol.

Only the standard library is required.
"""
from __future__ import annotations

import math
from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass
from typing import Any, ClassVar

from .geometry import BoundingBox
from .geometry_field import (
    DEFAULT_FRAME,
    DEFAULT_UNITS,
    IMPLICIT_FUNCTION,
    NEGATIVE_INSIDE,
    FieldError,
    GeometryField,
    Point,
    _points,
)

#: ``min`` of the operands: the points inside either one.
UNION = "union"
#: ``max`` of the operands: the points inside both.
INTERSECTION = "intersection"
#: ``max(left, -right)``: the points inside the left operand and outside the right.
DIFFERENCE = "difference"
OPERATORS = (UNION, INTERSECTION, DIFFERENCE)

#: Prefix separator used to address an operand's own controls through
#: :meth:`BooleanField.with_state` or :meth:`BooleanField.state` (``left__radius``).
STATE_SEPARATOR = "__"


def _operand(name: str, field: Any) -> GeometryField:
    if not isinstance(field, GeometryField):
        raise FieldError(f"{name} must be a GeometryField, got {type(field).__name__}")
    return field


def _check_operator(operator: str) -> str:
    if operator not in OPERATORS:
        raise FieldError(f"operator must be one of {OPERATORS}, got {operator!r}")
    return operator


@dataclass(frozen=True)
class BooleanField(GeometryField):
    """A CSG combination of two fields, as a field.

    Construct through :func:`union`, :func:`intersection` or :func:`difference`;
    the constructors are the API and they also handle the n-ary convenience form.
    The dataclass is frozen, so an edit returns a new field and the left and right
    operands are shared, not copied.
    """

    left: GeometryField
    right: GeometryField
    operator: str = UNION
    units: str = DEFAULT_UNITS
    coordinate_frame: str = DEFAULT_FRAME

    kind: ClassVar[str] = "boolean"
    quantity: ClassVar[str] = IMPLICIT_FUNCTION
    inside_sign: ClassVar[str] = NEGATIVE_INSIDE

    def __post_init__(self) -> None:
        object.__setattr__(self, "left", _operand("left", self.left))
        object.__setattr__(self, "right", _operand("right", self.right))
        object.__setattr__(self, "operator", _check_operator(self.operator))

    # -- evaluation -------------------------------------------------------- #

    def _operands(
        self, points: Iterable[Sequence[float]]
    ) -> tuple[list[Point], list[float], list[float]]:
        batch = _points(points)
        if not batch:
            return [], [], []
        return batch, self.left.signed_distance(batch), self.right.signed_distance(batch)

    def signed_distance(self, points: Iterable[Sequence[float]]) -> list[float]:
        _, left, right = self._operands(points)
        if self.operator == UNION:
            return [a if a <= b else b for a, b in zip(left, right, strict=True)]
        if self.operator == INTERSECTION:
            return [a if a >= b else b for a, b in zip(left, right, strict=True)]
        return [a if a >= -b else -b for a, b in zip(left, right, strict=True)]

    def gradient(self, points: Iterable[Sequence[float]]) -> list[Point]:
        """The gradient of the operand whose expression is active at each point.

        Away from the seam this is exact (the composition is that operand locally).
        On the seam the composed surface is not differentiable where it is not
        smooth, and one-sided derivatives are all a field can offer there; the
        caller that needs a seam normal should use the surface it extracted rather
        than this method. Ties go to the left operand, deterministically.
        """
        batch, left, right = self._operands(points)
        if not batch:
            return []
        left_gradient = self.left.gradient(batch)
        right_gradient = self.right.gradient(batch)
        result: list[Point] = []
        for index, (a, b) in enumerate(zip(left, right, strict=True)):
            if self.operator == UNION:
                use_left = a <= b
            elif self.operator == INTERSECTION:
                use_left = a >= b
            else:
                use_left = a >= -b
            if use_left:
                result.append(left_gradient[index])
            elif self.operator == DIFFERENCE:
                g = right_gradient[index]
                result.append((-g[0], -g[1], -g[2]))
            else:
                result.append(right_gradient[index])
        return result

    def bounds(self) -> BoundingBox | None:
        """Where the composed surface can be.

        ``union``: the box around both operands (``None`` if either is unbounded).
        ``intersection``: the overlap of the two boxes, or ``None`` if either operand
        is unbounded. Operands whose boxes do not overlap cannot intersect, and that
        is reported as a :class:`FieldError` rather than as a degenerate grid --
        sampling an empty sampling box would silently return "no surface".
        ``difference``: the left operand's box, which is a superset (a difference
        never grows the body).
        """
        left = self.left.bounds()
        right = self.right.bounds()
        if self.operator == UNION:
            if left is None or right is None:
                return None
            return BoundingBox(
                tuple(min(left.minimum[axis], right.minimum[axis]) for axis in range(3)),
                tuple(max(left.maximum[axis], right.maximum[axis]) for axis in range(3)),
            )
        if left is None:
            return right
        if right is None or self.operator == DIFFERENCE:
            return left
        minimum = tuple(max(left.minimum[axis], right.minimum[axis]) for axis in range(3))
        maximum = tuple(min(left.maximum[axis], right.maximum[axis]) for axis in range(3))
        if any(minimum[axis] > maximum[axis] for axis in range(3)):
            raise FieldError(
                f"the two operands of an intersection do not overlap: {left.minimum}.."
                f"{left.maximum} vs {right.minimum}..{right.maximum}"
            )
        return BoundingBox(minimum, maximum)

    # -- protocol bookkeeping ---------------------------------------------- #

    def state(self) -> dict[str, float]:
        """The operands' numeric controls, prefixed ``left__`` / ``right__``."""
        values: dict[str, float] = {}
        for name, operand in (("left", self.left), ("right", self.right)):
            for key, value in operand.state().items():
                values[f"{name}{STATE_SEPARATOR}{key}"] = value
        return values

    def with_state(self, **values: float) -> BooleanField:
        """Edit an operand's controls in place of the composition's own.

        A composition has no numeric controls of its own -- its state *is* its
        operands' state -- so the keys are addressed through the separator:
        ``field.with_state(left__radius=0.6)``. Unknown keys are an error, as they
        are for every other field, and the original is untouched.
        """
        known = set(self.state())
        unknown = set(values) - known
        if unknown:
            raise FieldError(
                f"{self.kind} has no control(s) {sorted(unknown)}; known controls are "
                f"{sorted(known) or 'none (the operands expose no numeric state)'}"
            )
        grouped: dict[str, dict[str, float]] = {"left": {}, "right": {}}
        for key, value in values.items():
            name, _, control = key.partition(STATE_SEPARATOR)
            grouped[name][control] = float(value)
        return BooleanField(
            self.left.with_state(**grouped["left"]) if grouped["left"] else self.left,
            self.right.with_state(**grouped["right"]) if grouped["right"] else self.right,
            self.operator,
            self.units,
            self.coordinate_frame,
        )

    def parameters(self) -> dict[str, Any]:
        return {
            "operator": self.operator,
            "left": self.left.parameters(),
            "right": self.right.parameters(),
        }

    def provenance(self) -> dict[str, Any]:
        return {
            "producer": "TensorCAD",
            "representation": "analytic",
            "kind": self.kind,
            "operator": self.operator,
        }

    def operands(self) -> tuple[GeometryField, GeometryField]:
        """The two operands, in the order ``signed_distance`` evaluates them."""
        return (self.left, self.right)

    def numpy_callable(self) -> Callable[[Any], Any] | None:
        """Vectorised composition, when both operands offer a vectorised form.

        ``minimum``/``maximum`` select an operand's value, they do not recompute it,
        so this path is bit-identical to the pure-Python one rather than close to it.
        """
        left = self.left.numpy_callable()
        right = self.right.numpy_callable()
        if left is None or right is None:
            return None
        operator = self.operator

        def evaluate(points: Any) -> Any:
            import numpy as np

            a = left(points)
            b = right(points)
            if operator == UNION:
                return np.minimum(a, b)
            if operator == INTERSECTION:
                return np.maximum(a, b)
            return np.maximum(a, -b)

        return evaluate


def _fold(operator: str, fields: Sequence[GeometryField], *, units: str,
          coordinate_frame: str) -> GeometryField:
    if not fields:
        raise FieldError(f"{operator} needs at least one field")
    result = _operand("field", fields[0])
    for index, field in enumerate(fields[1:], start=1):
        result = BooleanField(result, field, operator, units, coordinate_frame)
    return result


def union(*fields: GeometryField, units: str = DEFAULT_UNITS,
          coordinate_frame: str = DEFAULT_FRAME) -> GeometryField:
    """The solid inside **any** operand (``min`` of their fields).

    ``union(a)`` is ``a``; the n-ary form folds left, so
    ``union(a, b, c) == union(union(a, b), c)``.
    """
    return _fold(UNION, fields, units=units, coordinate_frame=coordinate_frame)


def intersection(*fields: GeometryField, units: str = DEFAULT_UNITS,
                 coordinate_frame: str = DEFAULT_FRAME) -> GeometryField:
    """The solid inside **every** operand (``max`` of their fields)."""
    return _fold(INTERSECTION, fields, units=units, coordinate_frame=coordinate_frame)


def difference(*fields: GeometryField, units: str = DEFAULT_UNITS,
               coordinate_frame: str = DEFAULT_FRAME) -> GeometryField:
    """The left operand with the others removed (``max(a, -b, -c, ...)``)."""
    return _fold(DIFFERENCE, fields, units=units, coordinate_frame=coordinate_frame)


def distance_to_surface(
    field: GeometryField, points: Iterable[Sequence[float]], *,
    rings: int = 64, meridians: int = 128,
) -> list[float]:
    """Brute-force unsigned distance to the surface, for validation only.

    Only meaningful for a composite whose *operands* have closed-form surfaces, so
    the sample cloud below assumes spheres: the point set is ``rings x meridians``
    points per sphere operand, and the minimum distance to that cloud is an upper
    bound on the true distance that tightens as the sampling densifies. The test
    suite uses it to check the conservativeness claim in the module docstring
    (``|value| >= distance``), not as a production query.
    """
    batch = _points(points)
    if not batch:
        return []
    cloud: list[Point] = []
    for operand in _leaves(field):
        radius = getattr(operand, "radius", None)
        center = getattr(operand, "center", None)
        if radius is None or center is None:
            raise FieldError("distance_to_surface only supports sphere operands")
        cx, cy, cz = center
        for ring in range(1, rings):
            theta = math.pi * ring / rings
            for meridian in range(meridians):
                phi = 2.0 * math.pi * meridian / meridians
                cloud.append(
                    (
                        cx + radius * math.sin(theta) * math.cos(phi),
                        cy + radius * math.sin(theta) * math.sin(phi),
                        cz + radius * math.cos(theta),
                    )
                )
        cloud.append((cx, cy, cz + radius))
        cloud.append((cx, cy, cz - radius))
    distances: list[float] = []
    for point in batch:
        best = math.inf
        for sample in cloud:
            dx = point[0] - sample[0]
            dy = point[1] - sample[1]
            dz = point[2] - sample[2]
            squared = dx * dx + dy * dy + dz * dz
            if squared < best:
                best = squared
        distances.append(math.sqrt(best))
    return distances


def _leaves(field: GeometryField) -> tuple[GeometryField, ...]:
    """The primitive fields a composition is built from, in order."""
    if isinstance(field, BooleanField):
        return _leaves(field.left) + _leaves(field.right)
    return (field,)