"""Continuous geometry fields: one query protocol for analytic and learned geometry.

A geometry in TensorCAD is a *field*, not a mesh: a signed quantity is defined
everywhere in space, and the surface is its zero level set. A `.ngo` object is one
implementation (a learned field). This module is the other half of the contract --
exact analytic primitives, and the protocol every field has to satisfy:

    value  = field.signed_distance(points)   # batch, (n, 3) -> (n,)
    normal = field.gradient(points)          # d(value)/dx, unnormalised is wrong
    inside = field.contains(points)          # value < surface_isovalue
    bounds = field.bounds()                  # where the surface can be
    state  = field.with_state(radius=0.6)    # an edit, as a new field

Why a protocol instead of a class hierarchy per shape: the same code has to run
over an analytic sphere (exact, microsecond, no dependencies) and a learned
ellipsoid (approximate, ~3 % of the radius, needs torch). Making them
interchangeable is what turns "we can fit a shape" into "we can edit, validate and
hand a shape to a solver" -- and it is what makes the analytic primitives usable
as *reference implementations* when judging a learned one.

Two semantics are recorded rather than assumed, because conflating them produces
confident nonsense:

* ``quantity`` -- ``"signed_distance"`` means ``|value|`` is the Euclidean
  distance to the surface (plane, sphere, capped cylinder). ``"implicit_function"``
  means only the zero set and the sign are meaningful (an ellipsoid: its distance
  has no closed form). A benchmark may compare *surfaces*, never raw values,
  across the two.
* ``inside_sign`` -- negative inside, the same convention as the NGO format, so a
  field from either source orients an extracted mesh the same way.

Only the standard library is required. Sampling a field into a triangle mesh is
the one operation that needs NumPy, and it is imported lazily inside
:func:`to_mesh`.
"""
from __future__ import annotations

import math
from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass
from typing import Any, ClassVar

from .geometry import BoundingBox, Mesh

Point = tuple[float, float, float]

SIGNED_DISTANCE = "signed_distance"
IMPLICIT_FUNCTION = "implicit_function"
#: A signed distance truncated to a band around the surface: ``|value|`` is the
#: distance only inside that band (``field.truncation_distance``), so magnitudes
#: beyond it are not comparable with an exact SDF. The learned NGO field is one.
TRUNCATED_SIGNED_DISTANCE = "truncated_signed_distance"
NEGATIVE_INSIDE = "negative"
DEFAULT_FRAME = "object_local"
DEFAULT_UNITS = "m"


class FieldError(ValueError):
    """Raised when a field cannot be constructed, queried or sampled."""


def _points(points: Iterable[Sequence[float]]) -> list[Point]:
    """Validate and normalise a batch of query points."""
    result: list[Point] = []
    for index, point in enumerate(points):
        values = tuple(float(value) for value in point)
        if len(values) != 3:
            raise FieldError(f"point {index} has {len(values)} coordinates, expected 3")
        for axis, value in enumerate(values):
            if not math.isfinite(value):
                raise FieldError(f"point {index} coordinate {axis} is not finite: {value!r}")
        result.append((values[0], values[1], values[2]))
    return result


def _normalise(vector: Point) -> Point:
    length = math.sqrt(vector[0] ** 2 + vector[1] ** 2 + vector[2] ** 2)
    if not math.isfinite(length) or length <= 0.0:
        raise FieldError(f"direction must be a finite non-zero vector, got {vector!r}")
    return (vector[0] / length, vector[1] / length, vector[2] / length)


def _positive(name: str, value: float) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0.0:
        raise FieldError(f"{name} must be a finite positive number, got {value!r}")
    return number


def _finite_point(name: str, value: Sequence[float]) -> Point:
    if len(value) != 3:
        raise FieldError(f"{name} needs exactly 3 coordinates, got {len(value)}")
    point = (float(value[0]), float(value[1]), float(value[2]))
    if not all(math.isfinite(axis) for axis in point):
        raise FieldError(f"{name} must be finite, got {value!r}")
    return point


class GeometryField:
    """The query protocol every geometry representation implements.

    Subclasses fill in :meth:`signed_distance`, :meth:`gradient`, :meth:`bounds`,
    :meth:`state`, :meth:`with_state` and :meth:`parameters`. Everything else --
    containment, metadata, numerical gradients, the batch contract -- is defined
    once, here, so an analytic primitive and a learned field cannot drift apart.
    """

    kind: ClassVar[str] = "abstract"
    quantity: ClassVar[str] = IMPLICIT_FUNCTION
    inside_sign: ClassVar[str] = NEGATIVE_INSIDE
    surface_isovalue: ClassVar[float] = 0.0
    units: str = DEFAULT_UNITS
    coordinate_frame: str = DEFAULT_FRAME

    # -- required by every field ------------------------------------------- #

    def signed_distance(self, points: Iterable[Sequence[float]]) -> list[float]:
        """Evaluate the field at ``points``; the surface is where this is zero."""
        raise NotImplementedError

    def gradient(self, points: Iterable[Sequence[float]]) -> list[Point]:
        """Analytic (or documented approximate) gradient, same order as ``points``."""
        raise NotImplementedError

    def bounds(self) -> BoundingBox | None:
        """Axis-aligned region that contains the surface, or ``None`` if unbounded."""
        return None

    def state(self) -> dict[str, float]:
        """The editable controls and their current values."""
        return {}

    def with_state(self, **values: float) -> GeometryField:
        """Return a new field with some controls changed; nothing is mutated."""
        raise NotImplementedError

    def parameters(self) -> dict[str, Any]:
        """The full construction parameter set, for the reproducibility record."""
        return {}

    # -- defined once, for every field ------------------------------------- #

    def __call__(self, points: Iterable[Sequence[float]]) -> list[float]:
        return self.signed_distance(points)

    def contains(self, points: Iterable[Sequence[float]]) -> list[bool]:
        """Inside/outside by the declared sign convention."""
        values = self.signed_distance(points)
        if self.inside_sign == NEGATIVE_INSIDE:
            return [value < self.surface_isovalue for value in values]
        return [value > self.surface_isovalue for value in values]

    def numerical_gradient(
        self, points: Iterable[Sequence[float]], *, step: float | None = None
    ) -> list[Point]:
        """Central differences, for fields that cannot write a gradient down.

        A learned field has no closed form and torch autograd is not always
        available (the extraction path is NumPy-only), so this is the fallback the
        protocol ships with. The step defaults to ``1e-6`` scaled by the field's
        own extent, which keeps truncation and round-off balanced for a field
        evaluated in ``float64``. A ``float32`` field must pass a coarser step --
        differencing two nearly equal float32 values below ~``1e-4`` of the
        object's size measures rounding, not a gradient.
        """
        batch = _points(points)
        if step is None:
            box = self.bounds()
            extent = max(box.sizes) if box is not None else 1.0
            step = max(1e-9, 1e-6 * extent)
        derivatives: list[list[float]] = []
        for axis in range(3):
            forward: list[Point] = []
            backward: list[Point] = []
            for point in batch:
                shifted = [point[0], point[1], point[2]]
                shifted[axis] += step
                forward.append((shifted[0], shifted[1], shifted[2]))
                shifted[axis] -= 2.0 * step
                backward.append((shifted[0], shifted[1], shifted[2]))
            plus = self.signed_distance(forward)
            minus = self.signed_distance(backward)
            derivatives.append(
                [(high - low) / (2.0 * step) for high, low in zip(plus, minus, strict=True)]
            )
        return [
            (derivatives[0][index], derivatives[1][index], derivatives[2][index])
            for index in range(len(batch))
        ]

    def metadata(self) -> dict[str, Any]:
        """Everything a consumer needs to reproduce and interpret this field."""
        box = self.bounds()
        return {
            "kind": self.kind,
            "quantity": self.quantity,
            "inside_sign": self.inside_sign,
            "surface_isovalue": self.surface_isovalue,
            "units": self.units,
            "coordinate_frame": self.coordinate_frame,
            "bounds": None if box is None else {key: list(value) for key, value in box.to_dict().items()},
            "state": self.state(),
            "parameters": self.parameters(),
            "provenance": dict(self.provenance()),
        }

    def provenance(self) -> dict[str, Any]:
        return {"producer": "TensorCAD", "representation": "analytic", "kind": self.kind}

    def numpy_callable(self) -> Callable[[Any], Any] | None:
        """Vectorised evaluator for grid sampling, or ``None`` when unavailable.

        Subclasses that can express their formula with NumPy return a callable
        taking an ``(n, 3)`` array. Sampling a grid point-by-point in Python is
        correct but ~50x slower, and an extraction is 10^5--10^6 queries, so the
        fast path matters; the pure-Python path stays the contract.
        """
        return None

    def __repr__(self) -> str:  # pragma: no cover - debugging aid
        return f"{type(self).__name__}({self.parameters()})"


@dataclass(frozen=True)
class PlaneField(GeometryField):
    """An oriented plane, ``n . x = offset``: the simplest exact reference field.

    Unbounded, so it has no bounds and cannot be sampled into a closed mesh -- it
    exists to pin the protocol's conventions (sign, gradient, state editing) to
    something whose answer is obvious.
    """

    normal: Point = (0.0, 0.0, 1.0)
    offset: float = 0.0

    kind: ClassVar[str] = "plane"
    quantity: ClassVar[str] = SIGNED_DISTANCE

    def __post_init__(self) -> None:
        object.__setattr__(self, "normal", _normalise(_finite_point("normal", self.normal)))
        if not math.isfinite(float(self.offset)):
            raise FieldError(f"offset must be finite, got {self.offset!r}")

    def signed_distance(self, points: Iterable[Sequence[float]]) -> list[float]:
        nx, ny, nz = self.normal
        return [nx * x + ny * y + nz * z - self.offset for x, y, z in _points(points)]

    def gradient(self, points: Iterable[Sequence[float]]) -> list[Point]:
        return [self.normal for _ in _points(points)]

    def state(self) -> dict[str, float]:
        return {"offset": self.offset}

    def with_state(self, **values: float) -> PlaneField:
        unknown = set(values) - {"offset"}
        if unknown:
            raise FieldError(f"{self.kind} has no control(s) {sorted(unknown)}")
        return PlaneField(self.normal, float(values.get("offset", self.offset)))

    def parameters(self) -> dict[str, Any]:
        return {"normal": list(self.normal), "offset": self.offset}

    def numpy_callable(self) -> Callable[[Any], Any] | None:
        try:
            import numpy as np
        except ImportError:  # pragma: no cover - numpy is in the test extra
            return None
        normal = np.asarray(self.normal, dtype=np.float64)
        offset = self.offset

        def evaluate(points: Any) -> Any:
            return np.asarray(points, dtype=np.float64) @ normal - offset

        return evaluate


@dataclass(frozen=True)
class SphereField(GeometryField):
    """A sphere: ``|x - c| - r``. Exact distance, exact gradient."""

    radius: float = 0.5
    center: Point = (0.0, 0.0, 0.0)

    kind: ClassVar[str] = "sphere"
    quantity: ClassVar[str] = SIGNED_DISTANCE

    def __post_init__(self) -> None:
        object.__setattr__(self, "radius", _positive("radius", self.radius))
        object.__setattr__(self, "center", _finite_point("center", self.center))

    def signed_distance(self, points: Iterable[Sequence[float]]) -> list[float]:
        cx, cy, cz = self.center
        radius = self.radius
        return [
            math.sqrt((x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2) - radius
            for x, y, z in _points(points)
        ]

    def gradient(self, points: Iterable[Sequence[float]]) -> list[Point]:
        cx, cy, cz = self.center
        gradients: list[Point] = []
        for x, y, z in _points(points):
            dx, dy, dz = x - cx, y - cy, z - cz
            length = math.sqrt(dx * dx + dy * dy + dz * dz)
            if length <= 1e-15:  # the centre: every direction is a subgradient
                gradients.append((0.0, 0.0, 0.0))
            else:
                gradients.append((dx / length, dy / length, dz / length))
        return gradients

    def bounds(self) -> BoundingBox:
        cx, cy, cz = self.center
        r = self.radius
        return BoundingBox((cx - r, cy - r, cz - r), (cx + r, cy + r, cz + r))

    def state(self) -> dict[str, float]:
        return {"radius": self.radius}

    def with_state(self, **values: float) -> SphereField:
        unknown = set(values) - {"radius"}
        if unknown:
            raise FieldError(f"{self.kind} has no control(s) {sorted(unknown)}")
        # The centre is a parameter, not a control: editing the radius must not
        # move the sphere, which is exactly the bug this test exists to catch.
        return SphereField(float(values.get("radius", self.radius)), self.center)

    def parameters(self) -> dict[str, Any]:
        return {"radius": self.radius, "center": list(self.center)}

    def numpy_callable(self) -> Callable[[Any], Any] | None:
        try:
            import numpy as np
        except ImportError:  # pragma: no cover
            return None
        center = np.asarray(self.center, dtype=np.float64)
        radius = self.radius

        def evaluate(points: Any) -> Any:
            delta = np.asarray(points, dtype=np.float64) - center
            return np.sqrt(np.einsum("ij,ij->i", delta, delta)) - radius

        return evaluate


@dataclass(frozen=True)
class CylinderField(GeometryField):
    """A capped circular cylinder, oriented along ``axis`` through ``center``."""

    radius: float = 0.3
    height: float = 1.0
    center: Point = (0.0, 0.0, 0.0)
    axis: Point = (0.0, 0.0, 1.0)

    kind: ClassVar[str] = "cylinder"
    quantity: ClassVar[str] = SIGNED_DISTANCE

    def __post_init__(self) -> None:
        object.__setattr__(self, "radius", _positive("radius", self.radius))
        object.__setattr__(self, "height", _positive("height", self.height))
        object.__setattr__(self, "center", _finite_point("center", self.center))
        object.__setattr__(self, "axis", _normalise(_finite_point("axis", self.axis)))

    def _components(self, point: Sequence[float]) -> tuple[float, float, float, Point]:
        """Return ``(radial_distance, axial_distance, axial_coordinate, radial_dir)``."""
        cx, cy, cz = self.center
        ax, ay, az = self.axis
        dx, dy, dz = point[0] - cx, point[1] - cy, point[2] - cz
        axial = dx * ax + dy * ay + dz * az
        rx, ry, rz = dx - axial * ax, dy - axial * ay, dz - axial * az
        radial_length = math.sqrt(rx * rx + ry * ry + rz * rz)
        if radial_length <= 1e-15:
            radial_direction: Point = (0.0, 0.0, 0.0)
        else:
            radial_direction = (rx / radial_length, ry / radial_length, rz / radial_length)
        return radial_length - self.radius, abs(axial) - 0.5 * self.height, axial, radial_direction

    def signed_distance(self, points: Iterable[Sequence[float]]) -> list[float]:
        values: list[float] = []
        for point in _points(points):
            radial, axial, _, _ = self._components(point)
            outside = math.hypot(max(radial, 0.0), max(axial, 0.0))
            values.append(outside + min(max(radial, axial), 0.0))
        return values

    def gradient(self, points: Iterable[Sequence[float]]) -> list[Point]:
        ax, ay, az = self.axis
        gradients: list[Point] = []
        for point in _points(points):
            radial, axial, axial_coordinate, radial_direction = self._components(point)
            axial_direction = (
                ax if axial_coordinate >= 0.0 else -ax,
                ay if axial_coordinate >= 0.0 else -ay,
                az if axial_coordinate >= 0.0 else -az,
            )
            if radial > 0.0 and axial > 0.0:
                # The corner region is the distance to the rim circle.
                gx = radial_direction[0] * radial + axial_direction[0] * axial
                gy = radial_direction[1] * radial + axial_direction[1] * axial
                gz = radial_direction[2] * radial + axial_direction[2] * axial
                length = math.sqrt(gx * gx + gy * gy + gz * gz)
                if length <= 1e-15:
                    gradients.append((0.0, 0.0, 0.0))
                else:
                    gradients.append((gx / length, gy / length, gz / length))
            elif radial >= axial:
                gradients.append(radial_direction)
            else:
                gradients.append(axial_direction)
        return gradients

    def bounds(self) -> BoundingBox:
        cx, cy, cz = self.center
        ax, ay, az = self.axis
        half = 0.5 * self.height
        r = self.radius
        # Extent of a capped cylinder along each world axis.
        ex = half * abs(ax) + r * math.sqrt(max(0.0, 1.0 - ax * ax))
        ey = half * abs(ay) + r * math.sqrt(max(0.0, 1.0 - ay * ay))
        ez = half * abs(az) + r * math.sqrt(max(0.0, 1.0 - az * az))
        return BoundingBox((cx - ex, cy - ey, cz - ez), (cx + ex, cy + ey, cz + ez))

    def state(self) -> dict[str, float]:
        return {"radius": self.radius, "height": self.height}

    def with_state(self, **values: float) -> CylinderField:
        unknown = set(values) - {"radius", "height"}
        if unknown:
            raise FieldError(f"{self.kind} has no control(s) {sorted(unknown)}")
        return CylinderField(
            float(values.get("radius", self.radius)),
            float(values.get("height", self.height)),
            self.center,
            self.axis,
        )

    def parameters(self) -> dict[str, Any]:
        return {
            "radius": self.radius,
            "height": self.height,
            "center": list(self.center),
            "axis": list(self.axis),
        }

    def numpy_callable(self) -> Callable[[Any], Any] | None:
        try:
            import numpy as np
        except ImportError:  # pragma: no cover
            return None
        center = np.asarray(self.center, dtype=np.float64)
        axis = np.asarray(self.axis, dtype=np.float64)
        radius = self.radius
        half = 0.5 * self.height

        def evaluate(points: Any) -> Any:
            delta = np.asarray(points, dtype=np.float64) - center
            axial = delta @ axis
            radial_vec = delta - axial[:, None] * axis
            radial = np.sqrt(np.einsum("ij,ij->i", radial_vec, radial_vec)) - radius
            along = np.abs(axial) - half
            outside = np.hypot(np.maximum(radial, 0.0), np.maximum(along, 0.0))
            return outside + np.minimum(np.maximum(radial, along), 0.0)

        return evaluate


@dataclass(frozen=True)
class EllipsoidField(GeometryField):
    """An axis-aligned ellipsoid, as an implicit function of the unit ball.

    ``f = sqrt(sum(((x - c)_i / a_i)^2)) - 1``. Its zero set and sign are exact;
    its magnitude is **not** a distance (an ellipsoid's distance has no closed
    form), which is why the field declares ``quantity = "implicit_function"``
    rather than pretending otherwise. Distance comparisons against a learned field
    therefore sample directions and measure surface position, never raw values.
    """

    semi_axes: Point = (0.5, 0.3, 0.2)
    center: Point = (0.0, 0.0, 0.0)

    kind: ClassVar[str] = "ellipsoid"
    quantity: ClassVar[str] = IMPLICIT_FUNCTION

    def __post_init__(self) -> None:
        axes = _finite_point("semi_axes", self.semi_axes)
        validated = tuple(_positive(f"semi_axes[{index}]", axis) for index, axis in enumerate(axes))
        object.__setattr__(self, "semi_axes", validated)
        object.__setattr__(self, "center", _finite_point("center", self.center))

    def _quadratic(self, point: Sequence[float]) -> tuple[float, Point]:
        cx, cy, cz = self.center
        ax, ay, az = self.semi_axes
        u = (point[0] - cx) / ax
        v = (point[1] - cy) / ay
        w = (point[2] - cz) / az
        return math.sqrt(u * u + v * v + w * w), (u, v, w)

    def signed_distance(self, points: Iterable[Sequence[float]]) -> list[float]:
        return [self._quadratic(point)[0] - 1.0 for point in _points(points)]

    def gradient(self, points: Iterable[Sequence[float]]) -> list[Point]:
        ax, ay, az = self.semi_axes
        gradients: list[Point] = []
        for point in _points(points):
            q, (u, v, w) = self._quadratic(point)
            if q <= 1e-15:
                gradients.append((0.0, 0.0, 0.0))
            else:
                gradients.append((u / (ax * q), v / (ay * q), w / (az * q)))
        return gradients

    def bounds(self) -> BoundingBox:
        cx, cy, cz = self.center
        ax, ay, az = self.semi_axes
        return BoundingBox((cx - ax, cy - ay, cz - az), (cx + ax, cy + ay, cz + az))

    def state(self) -> dict[str, float]:
        return {"semi_x": self.semi_axes[0], "semi_y": self.semi_axes[1], "semi_z": self.semi_axes[2]}

    def with_state(self, **values: float) -> EllipsoidField:
        unknown = set(values) - {"semi_x", "semi_y", "semi_z"}
        if unknown:
            raise FieldError(f"{self.kind} has no control(s) {sorted(unknown)}")
        current = self.state()
        current.update({key: float(value) for key, value in values.items()})
        return EllipsoidField(
            (current["semi_x"], current["semi_y"], current["semi_z"]), self.center
        )

    def parameters(self) -> dict[str, Any]:
        return {"semi_axes": list(self.semi_axes), "center": list(self.center)}

    def numpy_callable(self) -> Callable[[Any], Any] | None:
        try:
            import numpy as np
        except ImportError:  # pragma: no cover
            return None
        center = np.asarray(self.center, dtype=np.float64)
        axes = np.asarray(self.semi_axes, dtype=np.float64)

        def evaluate(points: Any) -> Any:
            q = (np.asarray(points, dtype=np.float64) - center) / axes
            return np.sqrt(np.einsum("ij,ij->i", q, q)) - 1.0

        return evaluate


def surface_point(field: GeometryField, direction: Sequence[float]) -> Point:
    """Walk from the origin-ish interior along ``direction`` to the surface.

    Used by the benchmark to compare *surface positions* between representations
    whose field magnitudes are not comparable (an exact SDF and a learned
    truncated field, or an implicit ellipsoid). The walk is a bisection on the
    sign change along the ray, so it needs only the protocol method every field
    has and inherits its accuracy from the field, not from extra assumptions.
    """
    box = field.bounds()
    if box is None:
        raise FieldError(f"{field.kind} is unbounded; a surface ray needs bounds")
    unit = _normalise(_finite_point("direction", direction))
    inside = box.minimum
    furthest = math.sqrt(
        sum(((max(abs(box.minimum[i]), abs(box.maximum[i]))) ** 2 for i in range(3)), 0.0)
    )
    origin = (0.0, 0.0, 0.0)
    if field.contains([origin])[0]:
        start = origin
    else:  # pragma: no cover - every primitive here contains its centre
        start = (
            0.5 * (inside[0] + box.maximum[0]),
            0.5 * (inside[1] + box.maximum[1]),
            0.5 * (inside[2] + box.maximum[2]),
        )
    low, high = 0.0, furthest
    for _ in range(80):
        middle = 0.5 * (low + high)
        probe = (start[0] + middle * unit[0], start[1] + middle * unit[1], start[2] + middle * unit[2])
        if field.contains([probe])[0]:
            low = middle
        else:
            high = middle
        if high - low <= 1e-12 * max(1.0, furthest):
            break
    radius = 0.5 * (low + high)
    return (start[0] + radius * unit[0], start[1] + radius * unit[1], start[2] + radius * unit[2])


@dataclass(frozen=True)
class SampledField:
    """The result of sampling a field into a mesh, with its evidence."""

    mesh: Mesh
    report: Any
    validation: Any
    quality: tuple[Any, ...]

    @property
    def passed(self) -> bool:
        return bool(self.validation.passed) and bool(self.report.watertight)

    def mesh_digest(self) -> str:
        """SHA-256 over the mesh's coordinates and indices, independent of order.

        Callers that hand a mesh to another tool (a solver, a report, a review)
        need to say *which* mesh they handed over; a file hash would do, but the
        mesh is in memory, so hash its exact float64 bits instead.
        """
        import hashlib
        import struct

        digest = hashlib.sha256()
        digest.update(struct.pack("<II", len(self.mesh.vertices), len(self.mesh.faces)))
        for vertex in self.mesh.vertices:
            digest.update(struct.pack("<3d", float(vertex[0]), float(vertex[1]), float(vertex[2])))
        for face in self.mesh.faces:
            digest.update(struct.pack("<3I", int(face[0]), int(face[1]), int(face[2])))
        return digest.hexdigest()

    def handoff(self) -> dict[str, Any]:
        """The record a downstream solver needs before it trusts this mesh.

        Deliberately self-describing: units, frame, digest, topology, measured
        volume and the quality numbers that are warnings rather than errors, so a
        finite-element code (TensorFEM) can accept or reject the surface on its own
        terms instead of re-deriving them.
        """
        return {
            "schema": "tensorcad.geometry.handoff/1",
            "mesh_sha256": self.mesh_digest(),
            "units": self.mesh.units,
            "coordinate_frame": self.mesh.coordinate_frame,
            "vertices": len(self.mesh.vertices),
            "triangles": len(self.mesh.faces),
            "watertight": bool(self.report.watertight),
            "edge_manifold": bool(self.report.edge_manifold),
            "euler_characteristic": int(self.report.euler_characteristic),
            "signed_volume_m3": float(self.report.signed_volume),
            "surface_area_m2": float(self.report.surface_area),
            "bounds_min_m": list(self.report.bounds_min),
            "bounds_max_m": list(self.report.bounds_max),
            "validation": dict(self.validation.to_dict()),
            "quality": {check.name: check.to_dict() for check in self.quality},
        }


#: Irrational fractions of one cell (1/phi**2, 1/phi**3, 1/phi**4) used to keep the
#: sampling grid off the surface on each axis independently.
_GRID_OFFSETS = (0.3819660112501051, 0.2360679774997897, 0.1458980337503154)


def to_mesh(
    field: GeometryField,
    resolution: int | Sequence[int] = 48,
    *,
    margin: int = 3,
    units: str | None = None,
    coordinate_frame: str | None = None,
    validate: bool = True,
) -> SampledField:
    """Sample ``field`` into a closed triangle mesh, then validate it.

    ``resolution`` is cells per axis, as one number or three. A slender body wants
    an anisotropic grid (a 8.6:1 hull at 64x64x64 wastes most of the cells and still
    under-resolves the length); the cell size is then whatever ``box / cells`` is on
    each axis.

    Sampling needs NumPy (imported here, not at module import), so the analytic
    fields stay usable in a standard-library-only environment as long as nobody
    asks for a mesh. The returned mesh goes through the same
    :func:`~tensorcad.validation.inspect_geometry` checks as the hull model --
    one acceptance pipeline for every representation, which is the whole point of
    the protocol.
    """
    box = field.bounds()
    if box is None:
        raise FieldError(f"{field.kind} is unbounded and cannot be sampled into a closed mesh")
    if min(box.sizes) <= 0.0:
        raise FieldError(f"{field.kind} has a degenerate bounds box: {box.sizes!r}")

    from .neural.extract import (  # optional NumPy path
        ExtractionError,
        cells_per_axis,
        extract_surface,
        summarise,
    )
    from .validation import inspect_geometry, inspect_mesh_quality

    evaluate = field.numpy_callable()

    def pure_python(points: Any) -> Any:  # pragma: no cover - every primitive ships NumPy
        return field.signed_distance(points)

    field_callable: Callable[[Any], Any] = evaluate if evaluate is not None else pure_python

    # Grid-vertex coincidence: a shape whose extent is a whole number of cells long
    # (r = 0.3 in a +-0.3 box, say) puts the level set exactly on grid vertices, and
    # marching tetrahedra then emits thousands of zero-area triangles -- 12.5 % of the
    # faces for a cylinder at 40 cells per axis. Shifting the sampling box by an
    # irrational fraction of one cell removes the coincidence without changing the
    # resolution, the units, or the margin. Each axis gets its own irrational offset.
    cells = cells_per_axis(resolution)
    pad = tuple(
        _GRID_OFFSETS[axis] * box.sizes[axis] / cells[axis] for axis in range(3)
    )
    minimum = tuple(box.minimum[axis] - pad[axis] for axis in range(3))
    maximum = tuple(box.maximum[axis] + pad[axis] for axis in range(3))

    try:
        mesh = extract_surface(
            field_callable,
            minimum,
            maximum,
            resolution,
            margin=margin,
            units=units or field.units,
            coordinate_frame=coordinate_frame or field.coordinate_frame,
        )
    except ExtractionError as exc:
        raise FieldError(f"sampling {field.kind} at {resolution}^3 produced no surface: {exc}") from exc

    report = summarise(mesh, resolution, expected_units=units or field.units)
    quality = tuple(inspect_mesh_quality(mesh))
    if validate:
        validation = _report(inspect_geometry(mesh) + list(quality))
    else:  # pragma: no cover - benchmark convenience
        from .validation import ValidationReport

        validation = ValidationReport(())
    return SampledField(mesh=mesh, report=report, validation=validation, quality=quality)


def _report(checks: list[Any]) -> Any:
    from .validation import ValidationReport

    return ValidationReport(tuple(checks))


__all__ = [
    "IMPLICIT_FUNCTION",
    "SIGNED_DISTANCE",
    "TRUNCATED_SIGNED_DISTANCE",
    "CylinderField",
    "EllipsoidField",
    "FieldError",
    "GeometryField",
    "PlaneField",
    "SampledField",
    "SphereField",
    "surface_point",
    "to_mesh",
]
