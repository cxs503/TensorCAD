"""STEP AP203/AP214 analytic B-rep reader: solids in, triangle mesh out.

Why this module exists
----------------------
:mod:`tensorcad.exchange` speaks ``FACETED_BREP``: a triangle soup that somebody
else already tessellated. Real CAD hands over ``ADVANCED_FACE`` /
``MANIFOLD_SOLID_BREP`` files instead, whose geometry is *analytic* -- planes,
cylinders and cones trimmed by ``EDGE_CURVE`` entities -- and those are the files a
shipyard, a class society or a CFD group actually attaches to a drawing. This
module reads that family and nothing more; it is the read side of the seam that
:mod:`tensorcad.adapters.ocp` covers on the write side when OCP is installed.

Method
------
1. parse the ISO-10303-21 exchange structure (comments, strings, complex
   instances, typed parameters) with no third-party parser;
2. resolve ``AXIS2_PLACEMENT_3D`` frames, points and surface/curve definitions;
3. walk each face's bounds: ``EDGE_LOOP`` -> ``ORIENTED_EDGE`` -> ``EDGE_CURVE``,
   sampling each edge between its two ``VERTEX_POINT``s in the direction the file
   declares (``ORIENTED_EDGE.orientation`` crossed with ``EDGE_CURVE.same_sense``)
   and mapping every sample into the surface's own ``(u, v)`` parameter space;
4. triangulate that parameter region -- ear clipping after hole bridging, so a
   face with holes is not "a fan over the concatenation" -- and subdivide each
   triangle in parameter space until its angular span is below
   ``angular_tolerance``, evaluating the exact surface at every new vertex;
5. weld by coordinate and measure: divergence-theorem volume, watertightness,
   Euler characteristic and a per-surface-type census.

What it will not do
-------------------
It will not guess. A surface type outside the supported table raises
:class:`StepUnsupportedStep` naming the entity, because a guessed tessellation is
indistinguishable from a correct one downstream and the reported volume would be
silently wrong. A face whose single closed boundary loop on a doubly curved
surface is ambiguous (two caps fit the same loop) is refused for the same reason.
Pass ``allow_unsupported=True`` for a partial mesh: the result is then marked
``complete=False`` and its volume is reported as ``None``, never as a number.

Units: STEP carries a ``LENGTH_UNIT``; SI prefixes and ``CONVERSION_BASED_UNIT``
(inch/foot) are honoured, the scale (file units per metre) and its provenance are
recorded, and nothing is assumed unless the file says nothing -- then the
millimetre convention is recorded *as an assumption*, not as a fact.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Sequence

from .geometry import Mesh
from .mesh import is_watertight, total_signed_volume

__all__ = [
    "Part21",
    "StepBrepImport",
    "StepError",
    "StepParseError",
    "StepTessellationError",
    "StepUnsupportedStep",
    "census",
    "read_step_brep",
    "STEP_BREP_SCHEMA",
]

STEP_BREP_SCHEMA = "tensorcad.step_brep/1"

Vec = tuple[float, float, float]


# --------------------------------------------------------------------------- #
# errors
# --------------------------------------------------------------------------- #
class StepError(RuntimeError):
    """Base class for every failure this reader reports."""


class StepParseError(StepError):
    """The exchange structure is not a well-formed ISO-10303-21 file."""


class StepTessellationError(StepError):
    """A face's trim could not be triangulated without inventing area."""


class StepUnsupportedStep(StepError):
    """The file contains entities this reader will not tessellate by guesswork."""

    def __init__(self, source: str, unsupported: dict[str, str]):
        self.unsupported = dict(unsupported)
        detail = "; ".join(f"{name} ({reason})" for name, reason in sorted(unsupported.items()))
        super().__init__(
            f"{source}: unsupported STEP entities: {detail}. "
            'An analytic read beyond this table needs OpenCASCADE (pip install "tensorcad[ocp]"); '
            "guessing here would audit the tessellator, not the model"
        )


# --------------------------------------------------------------------------- #
# ISO-10303-21 exchange structure
# --------------------------------------------------------------------------- #
class Ref(str):
    """An entity reference ``#12`` (the string holds the id without ``#``)."""


class Enum(str):
    """An enumeration value such as ``.T.`` or ``.MILLI.`` (stored upper-case)."""


@dataclass(frozen=True)
class Typed:
    """A typed parameter such as ``LENGTH_MEASURE(25.4)``."""

    keyword: str
    values: tuple[Any, ...]


def _strip_comments(text: str) -> str:
    out: list[str] = []
    i = 0
    count = len(text)
    in_string = False
    while i < count:
        ch = text[i]
        if in_string:
            out.append(ch)
            if ch == "'":
                if i + 1 < count and text[i + 1] == "'":
                    out.append("'")
                    i += 2
                    continue
                in_string = False
            i += 1
            continue
        if ch == "'":
            in_string = True
            out.append(ch)
            i += 1
            continue
        if ch == "/" and i + 1 < count and text[i + 1] == "*":
            end = text.find("*/", i + 2)
            if end < 0:
                break
            out.append(" ")
            i = end + 2
            continue
        out.append(ch)
        i += 1
    return "".join(out)


def _split_statements(text: str) -> list[str]:
    parts: list[str] = []
    depth = 0
    in_string = False
    current: list[str] = []
    i = 0
    count = len(text)
    while i < count:
        ch = text[i]
        if in_string:
            current.append(ch)
            if ch == "'":
                if i + 1 < count and text[i + 1] == "'":
                    current.append("'")
                    i += 2
                    continue
                in_string = False
            i += 1
            continue
        if ch == "'":
            in_string = True
            current.append(ch)
        elif ch == "(":
            depth += 1
            current.append(ch)
        elif ch == ")":
            depth -= 1
            current.append(ch)
        elif ch == ";" and depth == 0:
            statement = "".join(current).strip()
            if statement:
                parts.append(statement)
            current = []
        else:
            current.append(ch)
        i += 1
    tail = "".join(current).strip()
    if tail:
        parts.append(tail)
    return parts


def _split_top(text: str) -> list[str]:
    """Split a STEP argument list on top-level commas."""
    parts: list[str] = []
    depth = 0
    in_string = False
    current: list[str] = []
    i = 0
    count = len(text)
    while i < count:
        ch = text[i]
        if in_string:
            current.append(ch)
            if ch == "'":
                if i + 1 < count and text[i + 1] == "'":
                    current.append("'")
                    i += 2
                    continue
                in_string = False
            i += 1
            continue
        if ch == "'":
            in_string = True
            current.append(ch)
        elif ch == "(":
            depth += 1
            current.append(ch)
        elif ch == ")":
            depth -= 1
            current.append(ch)
        elif ch == "," and depth == 0:
            parts.append("".join(current).strip())
            current = []
        else:
            current.append(ch)
        i += 1
    last = "".join(current).strip()
    if last or parts:
        parts.append(last)
    return parts


_NUMBER = re.compile(r"^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$")
_TYPED = re.compile(r"^([A-Za-z_][A-Za-z_0-9]*)\s*\((.*)\)$", re.S)


def _parse_value(token: str) -> Any:
    text = token.strip()
    if not text:
        return None
    if text in ("$", "*"):
        return None
    if text.startswith("#"):
        return Ref(text[1:])
    if len(text) >= 2 and text.startswith("'") and text.endswith("'"):
        return text[1:-1].replace("''", "'")
    if len(text) > 2 and text.startswith(".") and text.endswith("."):
        return Enum(text.upper())
    if text.startswith("(") and text.endswith(")"):
        inner = text[1:-1]
        if not inner.strip():
            return []
        return [_parse_value(part) for part in _split_top(inner)]
    if _NUMBER.match(text):
        return float(text)
    match = _TYPED.match(text)
    if match is not None:
        return Typed(match.group(1).upper(), tuple(_parse_value(p) for p in _split_top(match.group(2))))
    return text


def _flag(value: Any) -> bool:
    """An optional STEP BOOLEAN: ``.T.`` true, ``.F.`` false, ``$`` (unset) true.

    ``Enum`` is a ``str`` subclass, so ``bool(Enum('.F.'))`` is ``True``: every
    boolean in this file must go through here, not through ``bool()``.
    """
    if value is None:
        return True
    return str(value).strip().upper() != ".F."


def _parse_value_list(raw: str) -> list[Any]:
    if not raw.strip():
        return []
    return [_parse_value(part) for part in _split_top(raw)]


_SIMPLE_STATEMENT = re.compile(r"^#(\d+)\s*=\s*([A-Za-z_][A-Za-z_0-9]*)\s*\((.*)\)$", re.S)
_COMPLEX_STATEMENT = re.compile(r"^#(\d+)\s*=\s*\((.*)\)$", re.S)


class Part21:
    """A parsed ISO-10303-21 data section: id -> (type, parsed arguments)."""

    def __init__(self, text: str):
        cleaned = _strip_comments(text)
        data_start = cleaned.find("DATA;")
        data_end = cleaned.rfind("ENDSEC;")
        if data_start < 0:
            raise StepParseError("no DATA; section: not an ISO-10303-21 exchange structure")
        body = cleaned[data_start + len("DATA;") : data_end if data_end > data_start else len(cleaned)]
        self.statements: dict[str, tuple[str, list[Any]]] = {}
        self.complex: dict[str, str] = {}
        for statement in _split_statements(body):
            if not statement.startswith("#"):
                continue
            match = _SIMPLE_STATEMENT.match(statement)
            if match is not None:
                entity_id, entity_type, raw = match.groups()
                self.statements[entity_id] = (entity_type.upper(), _parse_value_list(raw))
                continue
            complex_match = _COMPLEX_STATEMENT.match(statement)
            if complex_match is None:
                raise StepParseError(f"cannot parse STEP statement {statement[:80]!r}")
            entity_id, raw = complex_match.groups()
            self.complex[entity_id] = raw
            self.statements[entity_id] = ("", _parse_value_list(raw))
        if not self.statements:
            raise StepParseError("the DATA section contains no entity instances")
        self.header = _header_fields(cleaned)

    def type(self, entity_id: str) -> str:
        try:
            return self.statements[entity_id][0]
        except KeyError as error:
            raise StepParseError(f"reference #{entity_id} does not exist (truncated file?)") from error

    def args(self, entity_id: str) -> list[Any]:
        try:
            return self.statements[entity_id][1]
        except KeyError as error:
            raise StepParseError(f"reference #{entity_id} does not exist (truncated file?)") from error

    def of_type(self, *types: str) -> list[str]:
        wanted = set(types)
        return [entity_id for entity_id, (entity_type, _) in self.statements.items() if entity_type in wanted]

    def refs(self, value: Any) -> list[Ref]:
        if isinstance(value, Ref):
            return [value]
        if isinstance(value, list):
            out: list[Ref] = []
            for item in value:
                out.extend(self.refs(item))
            return out
        if isinstance(value, Typed):
            out = []
            for item in value.values:
                out.extend(self.refs(item))
            return out
        return []


def _header_fields(text: str) -> dict[str, Any]:
    fields: dict[str, Any] = {}
    schema = re.search(r"FILE_SCHEMA\s*\(\s*\((.*?)\)\s*\)", text, re.S)
    if schema is not None:
        fields["file_schema"] = " ".join(part.strip().strip("'") for part in _split_top(schema.group(1)))
    name = re.search(r"FILE_NAME\s*\(\s*(.*?)\)\s*;", text, re.S)
    if name is not None:
        parts = _split_top(name.group(1))
        if parts:
            fields["file_name"] = _parse_value(parts[0])
        if len(parts) >= 7:
            fields["preprocessor"] = _parse_value(parts[5])
            fields["originating_system"] = _parse_value(parts[6])
    return fields


# --------------------------------------------------------------------------- #
# small vector algebra
# --------------------------------------------------------------------------- #
def _add(a: Vec, b: Vec) -> Vec:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def _sub(a: Vec, b: Vec) -> Vec:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _mul(a: Vec, k: float) -> Vec:
    return (a[0] * k, a[1] * k, a[2] * k)


def _dot(a: Vec, b: Vec) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _cross(a: Vec, b: Vec) -> Vec:
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def _norm(a: Vec) -> float:
    return math.sqrt(_dot(a, a))


def _unit(a: Vec) -> Vec:
    length = _norm(a)
    if length == 0.0:
        raise StepParseError("zero-length direction in the STEP file")
    return (a[0] / length, a[1] / length, a[2] / length)


@dataclass(frozen=True)
class Frame:
    """A right-handed frame from an ``AXIS2_PLACEMENT_3D``."""

    origin: Vec
    z: Vec
    x: Vec
    y: Vec

    @staticmethod
    def build(origin: Vec, axis: Vec | None, ref_direction: Vec | None) -> "Frame":
        z = _unit(axis) if axis is not None else (0.0, 0.0, 1.0)
        helper = ref_direction if ref_direction is not None else (1.0, 0.0, 0.0)
        projected = _sub(helper, _mul(z, _dot(helper, z)))
        if _norm(projected) < 1e-9:
            # The file's ref_direction is parallel to the axis (DERIVED or degenerate).
            alternative = (1.0, 0.0, 0.0) if abs(z[0]) < 0.9 else (0.0, 1.0, 0.0)
            projected = _cross(z, alternative)
        x = _unit(projected)
        y = _cross(z, x)
        return Frame(origin=origin, z=z, x=x, y=y)


# --------------------------------------------------------------------------- #
# surfaces: evaluate / invert / periodicity / planarity
# --------------------------------------------------------------------------- #
class _Surface:
    kind = "SURFACE"
    periodic_u = False
    periodic_v = False
    planar = False
    #: True when the second parameter is an angle (sphere/torus) rather than a length
    v_angular = False

    def evaluate(self, u: float, v: float) -> Vec:  # pragma: no cover - interface
        raise NotImplementedError

    def inverse(self, point: Vec) -> tuple[float, float]:  # pragma: no cover - interface
        raise NotImplementedError

    def angular(self, du: float, dv: float) -> float:
        """Largest angular span of a parameter-space step, in radians."""
        return max(du if self.periodic_u else 0.0, dv if self.periodic_v else 0.0)


class Plane(_Surface):
    kind = "PLANE"
    planar = True

    def __init__(self, frame: Frame):
        self.frame = frame

    def evaluate(self, u: float, v: float) -> Vec:
        f = self.frame
        return _add(f.origin, _add(_mul(f.x, u), _mul(f.y, v)))

    def inverse(self, point: Vec) -> tuple[float, float]:
        f = self.frame
        d = _sub(point, f.origin)
        return (_dot(d, f.x), _dot(d, f.y))


class Cylinder(_Surface):
    kind = "CYLINDRICAL_SURFACE"
    periodic_u = True

    def __init__(self, frame: Frame, radius: float):
        self.frame = frame
        self.radius = radius

    def evaluate(self, u: float, v: float) -> Vec:
        f = self.frame
        radial = _add(_mul(f.x, math.cos(u)), _mul(f.y, math.sin(u)))
        return _add(f.origin, _add(_mul(radial, self.radius), _mul(f.z, v)))

    def inverse(self, point: Vec) -> tuple[float, float]:
        f = self.frame
        d = _sub(point, f.origin)
        angle = math.atan2(_dot(d, f.y), _dot(d, f.x)) % (2.0 * math.pi)
        return (angle, _dot(d, f.z))


class Cone(_Surface):
    kind = "CONICAL_SURFACE"
    periodic_u = True

    def __init__(self, frame: Frame, radius: float, semi_angle: float):
        self.frame = frame
        self.radius = radius
        self.semi_angle = semi_angle
        self.slope = math.tan(semi_angle)

    def radius_at(self, v: float) -> float:
        return self.radius + v * self.slope

    def evaluate(self, u: float, v: float) -> Vec:
        f = self.frame
        radial = _add(_mul(f.x, math.cos(u)), _mul(f.y, math.sin(u)))
        return _add(f.origin, _add(_mul(radial, self.radius_at(v)), _mul(f.z, v)))

    def inverse(self, point: Vec) -> tuple[float, float]:
        f = self.frame
        d = _sub(point, f.origin)
        angle = math.atan2(_dot(d, f.y), _dot(d, f.x)) % (2.0 * math.pi)
        return (angle, _dot(d, f.z))


class Sphere(_Surface):
    kind = "SPHERICAL_SURFACE"
    periodic_u = True
    v_angular = True

    def __init__(self, frame: Frame, radius: float):
        self.frame = frame
        self.radius = radius

    def evaluate(self, u: float, v: float) -> Vec:
        f = self.frame
        cos_v = math.cos(v)
        radial = _add(_mul(f.x, cos_v * math.cos(u)), _mul(f.y, cos_v * math.sin(u)))
        return _add(f.origin, _add(_mul(radial, self.radius), _mul(f.z, self.radius * math.sin(v))))

    def inverse(self, point: Vec) -> tuple[float, float]:
        f = self.frame
        d = _sub(point, f.origin)
        angle = math.atan2(_dot(d, f.y), _dot(d, f.x)) % (2.0 * math.pi)
        height = max(-1.0, min(1.0, _dot(d, f.z) / self.radius))
        return (angle, math.asin(height))


class Torus(_Surface):
    kind = "TOROIDAL_SURFACE"
    periodic_u = True
    periodic_v = True
    v_angular = True

    def __init__(self, frame: Frame, major: float, minor: float):
        self.frame = frame
        self.major = major
        self.minor = minor

    def evaluate(self, u: float, v: float) -> Vec:
        f = self.frame
        radial = _add(_mul(f.x, math.cos(u)), _mul(f.y, math.sin(u)))
        tube = self.major + self.minor * math.cos(v)
        return _add(f.origin, _add(_mul(radial, tube), _mul(f.z, self.minor * math.sin(v))))

    def inverse(self, point: Vec) -> tuple[float, float]:
        f = self.frame
        d = _sub(point, f.origin)
        angle = math.atan2(_dot(d, f.y), _dot(d, f.x)) % (2.0 * math.pi)
        planar = math.hypot(_dot(d, f.x), _dot(d, f.y))
        return (angle, math.atan2(_dot(d, f.z), planar - self.major))


# --------------------------------------------------------------------------- #
# curves: used to sample an edge between its two vertices
# --------------------------------------------------------------------------- #
class _Curve:
    periodic = False
    domain: tuple[float, float] = (0.0, 0.0)

    def evaluate(self, t: float) -> Vec:  # pragma: no cover - interface
        raise NotImplementedError

    def inverse(self, point: Vec) -> float:  # pragma: no cover - interface
        raise NotImplementedError


class LineCurve(_Curve):
    def __init__(self, origin: Vec, direction: Vec):
        self.origin = origin
        self.direction = _unit(direction)

    def evaluate(self, t: float) -> Vec:
        return _add(self.origin, _mul(self.direction, t))

    def inverse(self, point: Vec) -> float:
        return _dot(_sub(point, self.origin), self.direction)


class CircleCurve(_Curve):
    periodic = True

    def __init__(self, frame: Frame, radius: float):
        self.frame = frame
        self.radius = radius

    def evaluate(self, t: float) -> Vec:
        f = self.frame
        return _add(
            f.origin,
            _add(_mul(f.x, self.radius * math.cos(t)), _mul(f.y, self.radius * math.sin(t))),
        )

    def inverse(self, point: Vec) -> float:
        f = self.frame
        d = _sub(point, f.origin)
        return math.atan2(_dot(d, f.y), _dot(d, f.x)) % (2.0 * math.pi)


class EllipseCurve(_Curve):
    periodic = True

    def __init__(self, frame: Frame, semi_x: float, semi_y: float):
        self.frame = frame
        self.semi_x = semi_x
        self.semi_y = semi_y

    def evaluate(self, t: float) -> Vec:
        f = self.frame
        return _add(
            f.origin,
            _add(_mul(f.x, self.semi_x * math.cos(t)), _mul(f.y, self.semi_y * math.sin(t))),
        )

    def inverse(self, point: Vec) -> float:
        f = self.frame
        d = _sub(point, f.origin)
        return math.atan2(_dot(d, f.y) / self.semi_y, _dot(d, f.x) / self.semi_x) % (2.0 * math.pi)


class BSplineCurve(_Curve):
    """Non-rational or rational B-spline curve, evaluated by de Boor."""

    def __init__(self, degree: int, points: Sequence[Vec], knots: Sequence[float], weights: Sequence[float] | None):
        self.degree = int(degree)
        self.points = list(points)
        self.knots = [float(k) for k in knots]
        self.weights = [1.0] * len(self.points) if weights is None else [float(w) for w in weights]
        if len(self.knots) != len(self.points) + self.degree + 1:
            raise StepParseError(
                f"B-spline curve has {len(self.points)} control points, degree {self.degree} "
                f"and {len(self.knots)} knots -- expected {len(self.points) + self.degree + 1}"
            )
        self.domain = (self.knots[self.degree], self.knots[len(self.points)])

    def evaluate(self, t: float) -> Vec:
        p = self.degree
        low, high = self.domain
        t = min(max(t, low), high)
        count = len(self.points)
        span = p
        while span < count - 1 and t >= self.knots[span + 1]:
            span += 1
        d = [
            (
                self.points[j][0] * self.weights[j],
                self.points[j][1] * self.weights[j],
                self.points[j][2] * self.weights[j],
                self.weights[j],
            )
            for j in range(span - p, span + 1)
        ]
        for r in range(1, p + 1):
            for j in range(p, r - 1, -1):
                i = span - p + j
                denominator = self.knots[i + p - r + 1] - self.knots[i]
                alpha = 0.0 if denominator == 0.0 else (t - self.knots[i]) / denominator
                previous = d[j - 1]
                current = d[j]
                d[j] = tuple((1.0 - alpha) * previous[k] + alpha * current[k] for k in range(4))
        x, y, z, w = d[p]
        if w == 0.0:
            raise StepParseError("B-spline evaluation hit a zero weight")
        return (x / w, y / w, z / w)

    def inverse(self, point: Vec) -> float:
        low, high = self.domain
        samples = 400
        best_t = low
        best = float("inf")
        for index in range(samples + 1):
            t = low + (high - low) * index / samples
            distance = _norm(_sub(self.evaluate(t), point))
            if distance < best:
                best = distance
                best_t = t
        step = (high - low) / samples
        left = max(low, best_t - step)
        right = min(high, best_t + step)
        for _ in range(60):
            m1 = left + (right - left) / 3.0
            m2 = right - (right - left) / 3.0
            if _norm(_sub(self.evaluate(m1), point)) < _norm(_sub(self.evaluate(m2), point)):
                right = m2
            else:
                left = m1
        return 0.5 * (left + right)


# --------------------------------------------------------------------------- #
# entity model: resolve references into the geometry above
# --------------------------------------------------------------------------- #
_SURFACE_REASONS = {
    "B_SPLINE_SURFACE_WITH_KNOTS": "B-spline surface inversion needs Newton iteration; out of scope for this reader",
    "B_SPLINE_SURFACE": "B-spline surface inversion needs Newton iteration; out of scope for this reader",
    "RATIONAL_B_SPLINE_SURFACE": "rational B-spline surface inversion is out of scope for this reader",
    "SURFACE_OF_REVOLUTION": "swept-surface inversion is out of scope for this reader",
    "SURFACE_OF_LINEAR_EXTRUSION": "swept-surface inversion is out of scope for this reader",
    "RECTANGULAR_TRIMMED_SURFACE": "composite trimmed surfaces are out of scope for this reader",
    "CURVE_BOUNDED_SURFACE": "composite trimmed surfaces are out of scope for this reader",
    "OFFSET_SURFACE": "offset-surface inversion is out of scope for this reader",
}

_CURVE_REASONS = {
    "POLYLINE": "POLYLINE edges are not sampled; write a B_SPLINE_CURVE or a LINE instead",
    "TRIMMED_CURVE": "TRIMMED_CURVE edges are not unwound; write the basis curve instead",
    "COMPOSITE_CURVE": "COMPOSITE_CURVE edges are not unwound",
    "SURFACE_CURVE": "SURFACE_CURVE edges are not unwound; use the 3D curve",
    "SEAM_CURVE": "SEAM_CURVE edges are not unwound; use the 3D curve",
}

SUPPORTED_SURFACES = (
    "PLANE",
    "CYLINDRICAL_SURFACE",
    "CONICAL_SURFACE",
    "SPHERICAL_SURFACE",
    "TOROIDAL_SURFACE",
)
SUPPORTED_CURVES = ("LINE", "CIRCLE", "ELLIPSE", "B_SPLINE_CURVE_WITH_KNOTS")


class _Model:
    def __init__(self, part: Part21, scale: float):
        self.part = part
        self.scale = scale  # file units per metre
        self._points: dict[str, Vec] = {}
        self._directions: dict[str, Vec] = {}
        self._frames: dict[str, Frame] = {}
        self._surfaces: dict[str, _Surface] = {}
        self._curves: dict[str, _Curve] = {}

    def point(self, entity_id: str) -> Vec:
        cached = self._points.get(entity_id)
        if cached is not None:
            return cached
        entity_type = self.part.type(entity_id)
        if entity_type not in ("CARTESIAN_POINT", "POINT"):
            raise StepParseError(f"#{entity_id} is {entity_type}, not a CARTESIAN_POINT")
        coordinates = self.part.args(entity_id)[1]
        values = [float(value) for value in coordinates if isinstance(value, float)]
        while len(values) < 3:
            values.append(0.0)
        result = (values[0] / self.scale, values[1] / self.scale, values[2] / self.scale)
        self._points[entity_id] = result
        return result

    def direction(self, entity_id: str) -> Vec:
        cached = self._directions.get(entity_id)
        if cached is not None:
            return cached
        entity_type = self.part.type(entity_id)
        if entity_type not in ("DIRECTION", "VECTOR"):
            raise StepParseError(f"#{entity_id} is {entity_type}, not a DIRECTION")
        values = [float(value) for value in self.part.args(entity_id)[1] if isinstance(value, float)]
        while len(values) < 3:
            values.append(0.0)
        result = (values[0], values[1], values[2])
        self._directions[entity_id] = result
        return result

    def frame(self, entity_id: str) -> Frame:
        cached = self._frames.get(entity_id)
        if cached is not None:
            return cached
        entity_type = self.part.type(entity_id)
        if entity_type == "AXIS2_PLACEMENT_2D":
            args = self.part.args(entity_id)
            origin = self.point(self.part.refs(args[1])[0])
            ref = self.direction(self.part.refs(args[2])[0]) if args[2] else None
            result = Frame.build(origin, (0.0, 0.0, 1.0), ref)
        elif entity_type == "AXIS2_PLACEMENT_3D":
            args = self.part.args(entity_id)
            origin = self.point(self.part.refs(args[1])[0]) if args[1] else (0.0, 0.0, 0.0)
            axis = self.direction(self.part.refs(args[2])[0]) if len(args) > 2 and args[2] else None
            ref_direction = self.direction(self.part.refs(args[3])[0]) if len(args) > 3 and args[3] else None
            result = Frame.build(origin, axis, ref_direction)
        else:
            raise StepParseError(f"#{entity_id} is {entity_type}, not an AXIS2_PLACEMENT_3D")
        self._frames[entity_id] = result
        return result

    def surface(self, entity_id: str) -> _Surface:
        cached = self._surfaces.get(entity_id)
        if cached is not None:
            return cached
        entity_type = self.part.type(entity_id)
        args = self.part.args(entity_id)
        if entity_type == "PLANE":
            result: _Surface = Plane(self.frame(self.part.refs(args[1])[0]))
        elif entity_type == "CYLINDRICAL_SURFACE":
            result = Cylinder(self.frame(self.part.refs(args[1])[0]), float(args[2]) / self.scale)
        elif entity_type == "CONICAL_SURFACE":
            result = Cone(self.frame(self.part.refs(args[1])[0]), float(args[2]) / self.scale, float(args[3]))
        elif entity_type == "SPHERICAL_SURFACE":
            result = Sphere(self.frame(self.part.refs(args[1])[0]), float(args[2]) / self.scale)
        elif entity_type == "TOROIDAL_SURFACE":
            result = Torus(
                self.frame(self.part.refs(args[1])[0]),
                float(args[2]) / self.scale,
                float(args[3]) / self.scale,
            )
        else:
            raise StepUnsupportedStep(
                str(entity_id), {entity_type: _SURFACE_REASONS.get(entity_type, "surface type not implemented")}
            )
        self._surfaces[entity_id] = result
        return result

    def curve(self, entity_id: str) -> _Curve:
        cached = self._curves.get(entity_id)
        if cached is not None:
            return cached
        entity_type = self.part.type(entity_id)
        args = self.part.args(entity_id)
        if entity_type == "LINE":
            origin = self.point(self.part.refs(args[1])[0])
            vector_id = self.part.refs(args[2])[0]
            vector_args = self.part.args(vector_id)
            direction_ref = self.part.refs(vector_args)[0]
            magnitude = next((float(value) for value in vector_args if isinstance(value, float)), 1.0)
            direction = _mul(self.direction(direction_ref), magnitude)
            result: _Curve = LineCurve(origin, direction)
        elif entity_type == "CIRCLE":
            result = CircleCurve(self.frame(self.part.refs(args[1])[0]), float(args[2]) / self.scale)
        elif entity_type == "ELLIPSE":
            result = EllipseCurve(
                self.frame(self.part.refs(args[1])[0]),
                float(args[2]) / self.scale,
                float(args[3]) / self.scale,
            )
        elif entity_type in ("B_SPLINE_CURVE_WITH_KNOTS", "B_SPLINE_CURVE"):
            result = self._bspline_curve(args)
        elif entity_type in _CURVE_REASONS:
            raise StepUnsupportedStep(str(entity_id), {entity_type: _CURVE_REASONS[entity_type]})
        else:
            raise StepUnsupportedStep(str(entity_id), {entity_type: "curve type not implemented"})
        self._curves[entity_id] = result
        return result

    def _bspline_curve(self, args: list[Any]) -> _Curve:
        degree = int(args[1])
        points = [self.point(ref) for ref in self.part.refs(args[2])]
        weights: list[float] | None = None
        if len(args) > 8 and isinstance(args[8], list) and args[8]:
            weights = [float(value) for value in args[8]]
        multiplicities = [int(value) for value in (args[6] or [])]
        knots = [float(value) for value in (args[7] or [])]
        expanded: list[float] = []
        for knot, multiplicity in zip(knots, multiplicities):
            expanded.extend([knot / self.scale] * multiplicity)
        return BSplineCurve(degree, points, expanded, weights)

    def vertex_point(self, entity_id: str) -> Vec:
        args = self.part.args(entity_id)
        return self.point(self.part.refs(args[1])[0])# --------------------------------------------------------------------------- #
# triangulation in parameter space
# --------------------------------------------------------------------------- #
def _signed_area(polygon: Sequence[tuple[float, float]]) -> float:
    total = 0.0
    count = len(polygon)
    for index in range(count):
        x1, y1 = polygon[index]
        x2, y2 = polygon[(index + 1) % count]
        total += x1 * y2 - x2 * y1
    return 0.5 * total


def _cross2(a: tuple[float, float], b: tuple[float, float], c: tuple[float, float]) -> float:
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def _inside(point: tuple[float, float], a: tuple[float, float], b: tuple[float, float], c: tuple[float, float]) -> bool:
    d1 = _cross2(a, b, point)
    d2 = _cross2(b, c, point)
    d3 = _cross2(c, a, point)
    has_negative = d1 < -1e-14 or d2 < -1e-14 or d3 < -1e-14
    has_positive = d1 > 1e-14 or d2 > 1e-14 or d3 > 1e-14
    return not (has_negative and has_positive)


def _segments_intersect(
    a: tuple[float, float],
    b: tuple[float, float],
    c: tuple[float, float],
    d: tuple[float, float],
) -> bool:
    d1 = _cross2(a, b, c)
    d2 = _cross2(a, b, d)
    d3 = _cross2(c, d, a)
    d4 = _cross2(c, d, b)
    return (
        ((d1 > 1e-14 and d2 < -1e-14) or (d1 < -1e-14 and d2 > 1e-14))
        and ((d3 > 1e-14 and d4 < -1e-14) or (d3 < -1e-14 and d4 > 1e-14))
    )


def _visible(
    outer: Sequence[tuple[float, float]],
    hole: Sequence[tuple[float, float]],
    outer_index: int,
    hole_index: int,
) -> bool:
    """Is the bridge from a hole vertex to an outer vertex free of crossings?"""
    a = hole[hole_index]
    b = outer[outer_index]
    for rings, skip in ((outer, outer_index), (hole, hole_index)):
        count = len(rings)
        for index in range(count):
            nxt = (index + 1) % count
            if index == skip or nxt == skip:
                continue
            if _segments_intersect(a, b, rings[index], rings[nxt]):
                return False
    return True


def _bridge_holes(
    outer: list[tuple[float, float]], holes: Sequence[list[tuple[float, float]]]
) -> list[tuple[float, float]]:
    """Merge every hole into the outer loop through a visibility-checked bridge."""
    merged = list(outer)
    for hole in holes:
        if len(hole) < 3:
            continue
        start = max(range(len(hole)), key=lambda index: (hole[index][0], hole[index][1]))
        candidates = sorted(range(len(merged)), key=lambda index: (-merged[index][0], -abs(merged[index][1])))
        outer_index = None
        for candidate in candidates:
            if _visible(merged, hole, candidate, start):
                outer_index = candidate
                break
        if outer_index is None:
            outer_index = min(
                range(len(merged)),
                key=lambda index: (merged[index][0] - hole[start][0]) ** 2
                + (merged[index][1] - hole[start][1]) ** 2,
            )
        merged = merged[: outer_index + 1] + hole[start:] + hole[: start + 1] + merged[outer_index:]
    return merged


def _same_point(a: tuple[float, float], b: tuple[float, float], tolerance: float = 1e-15) -> bool:
    return abs(a[0] - b[0]) <= tolerance and abs(a[1] - b[1]) <= tolerance


def _ear_clip(polygon: Sequence[tuple[float, float]]) -> list[tuple[tuple[float, float], ...]] | None:
    """Triangulate a simple CCW polygon, or return ``None`` if no valid ear exists.

    Two details matter for a polygon produced by hole bridging. First, the bridge
    vertices appear twice with the same coordinates, so the "does this ear contain
    another vertex" test has to compare coordinates as well as indices -- otherwise
    every ear touching the bridge is rejected. Second, there is no "clip something to
    make progress" fallback: silently emitting an overlapping triangle is how a
    tessellator reports a volume that is 2.5x too large and still looks finished.
    ``None`` propagates to an honest :class:`StepTessellationError`.
    """
    indices = list(range(len(polygon)))
    triangles: list[tuple[tuple[float, float], ...]] = []
    guard = 0
    limit = 8 * len(indices) + 100
    while len(indices) > 3 and guard < limit:
        guard += 1
        clipped = False
        for position in range(len(indices)):
            i_prev = indices[position - 1]
            i_cur = indices[position]
            i_next = indices[(position + 1) % len(indices)]
            a, b, c = polygon[i_prev], polygon[i_cur], polygon[i_next]
            if _cross2(a, b, c) <= 1e-16:
                continue
            blocked = False
            for other in indices:
                if other in (i_prev, i_cur, i_next):
                    continue
                candidate = polygon[other]
                if _same_point(candidate, a) or _same_point(candidate, b) or _same_point(candidate, c):
                    continue
                if _inside(candidate, a, b, c):
                    blocked = True
                    break
            if blocked:
                continue
            triangles.append((a, b, c))
            del indices[position]
            clipped = True
            break
        if not clipped:
            return None
    if len(indices) == 3:
        triangles.append((polygon[indices[0]], polygon[indices[1]], polygon[indices[2]]))
    return triangles


def _dedupe(points: Sequence[tuple[float, float]], tolerance: float = 1e-12) -> list[tuple[float, float]]:
    out: list[tuple[float, float]] = []
    for point in points:
        if out and abs(point[0] - out[-1][0]) <= tolerance and abs(point[1] - out[-1][1]) <= tolerance:
            continue
        out.append((float(point[0]), float(point[1])))
    while len(out) > 1 and abs(out[0][0] - out[-1][0]) <= tolerance and abs(out[0][1] - out[-1][1]) <= tolerance:
        out.pop()
    return out


def _unwrap(points: Sequence[tuple[float, float]], surface: _Surface) -> list[tuple[float, float]]:
    """Make a parameter trace continuous across a periodic seam."""
    fixed: list[tuple[float, float]] = [(float(points[0][0]), float(points[0][1]))]
    for raw_u, raw_v in points[1:]:
        previous_u, previous_v = fixed[-1]
        u, v = float(raw_u), float(raw_v)
        if surface.periodic_u:
            while u - previous_u > math.pi:
                u -= 2.0 * math.pi
            while previous_u - u > math.pi:
                u += 2.0 * math.pi
        if surface.periodic_v:
            while v - previous_v > math.pi:
                v -= 2.0 * math.pi
            while previous_v - v > math.pi:
                v += 2.0 * math.pi
        fixed.append((u, v))
    return fixed


@dataclass
class _FaceResult:
    triangles: list[tuple[Vec, Vec, Vec]] = field(default_factory=list)
    surface_kind: str = ""
    depth_exhausted: int = 0
    ambiguous: str | None = None


def _refine(
    triangle: tuple[tuple[float, float], tuple[float, float], tuple[float, float]],
    surface: _Surface,
    tolerance: float,
    max_depth: int,
) -> tuple[list[tuple[tuple[float, float], tuple[float, float], tuple[float, float]]], int]:
    """Bisect a parameter triangle until no edge spans more than ``tolerance`` radians.

    Longest-edge bisection rather than 1-to-4 subdivision: a sliver that spans a whole
    revolution but only a millimetre of height needs ~2*pi/tolerance splits, not
    (2*pi/tolerance)^2, and a fan triangulation of a revolved trim produces exactly
    those slivers.
    """
    out: list[tuple[tuple[float, float], tuple[float, float], tuple[float, float]]] = []
    exhausted = 0
    stack: list[tuple[tuple[tuple[float, float], tuple[float, float], tuple[float, float]], int]] = [(triangle, 0)]
    while stack:
        current, depth = stack.pop()
        a, b, c = current
        spans = (
            surface.angular(abs(a[0] - b[0]), abs(a[1] - b[1])),
            surface.angular(abs(b[0] - c[0]), abs(b[1] - c[1])),
            surface.angular(abs(c[0] - a[0]), abs(c[1] - a[1])),
        )
        longest = max(spans)
        if surface.planar or longest <= tolerance or depth >= max_depth:
            if not surface.planar and longest > tolerance:
                exhausted += 1
            out.append(current)
            continue
        if spans[0] == longest:
            first, second = a, b
            remaining = c
        elif spans[1] == longest:
            first, second = b, c
            remaining = a
        else:
            first, second = c, a
            remaining = b
        middle = (0.5 * (first[0] + second[0]), 0.5 * (first[1] + second[1]))
        stack.append(((first, middle, remaining), depth + 1))
        stack.append(((middle, second, remaining), depth + 1))
    return out, exhausted


def _v_range_at(polygon: Sequence[tuple[float, float]], u: float) -> tuple[float, float] | None:
    """The v interval the parameter polygon covers at ``u`` (a vertical line query)."""
    values: list[float] = []
    count = len(polygon)
    for index in range(count):
        p = polygon[index]
        q = polygon[(index + 1) % count]
        du = q[0] - p[0]
        if abs(du) <= 1e-15:
            if abs(p[0] - u) <= 1e-12:
                values.extend((p[1], q[1]))
            continue
        if (p[0] - u) * (q[0] - u) <= 0.0:
            values.append(p[1] + (q[1] - p[1]) * (u - p[0]) / du)
    if len(values) < 2:
        return None
    low, high = min(values), max(values)
    if high - low <= 1e-15:
        return None
    return (low, high)


def _wrap_strip_triangles(
    surface: _Surface,
    polygon: Sequence[tuple[float, float]],
    tolerance: float,
) -> list[tuple[Vec, Vec, Vec]] | None:
    """Tessellate a trim that wraps a periodic surface exactly once as a strip grid.

    A full-revolution face is the common case for any shaft, bore or pipe, and its
    parameter region is a rectangle (or a near-uniform strip). Sampling that rectangle
    on a ``u`` grid is both exact along the generator and far cheaper than ear-clipping
    a fan of slivers and subdividing it. ``None`` means "not this shape", and the
    caller falls through to the general ear-clipping path.
    """
    if not surface.periodic_u:
        return None
    u_values = [point[0] for point in polygon]
    v_values = [point[1] for point in polygon]
    u_low, u_high = min(u_values), max(u_values)
    v_low, v_high = min(v_values), max(v_values)
    if u_low > 1e-9 or u_high < 2.0 * math.pi - 1e-9:
        return None
    width = v_high - v_low
    if width <= 0.0:
        return None
    area = abs(_signed_area(polygon))
    if abs(area - 2.0 * math.pi * width) > 0.01 * 2.0 * math.pi * width:
        return None
    # The u grid is the boundary's own sample set (not an independent uniform grid):
    # the rim vertices of this face must coincide with the vertices the neighbouring
    # planar/annular faces already put on the same edge, or the mesh has a crack.
    u_grid = sorted(set(round(point[0], 12) for point in polygon))
    if 0.0 not in u_grid:
        u_grid.insert(0, 0.0)
    if abs(u_grid[-1] - 2.0 * math.pi) > 1e-12:
        u_grid.append(2.0 * math.pi)
    if len(u_grid) < 9:
        steps = min(720, max(8, int(math.ceil(2.0 * math.pi / max(tolerance, 1e-6)))))
        u_grid = [2.0 * math.pi * index / steps for index in range(steps + 1)]
    steps_v = 1
    if surface.v_angular:
        steps_v = min(180, max(1, int(math.ceil(width / max(tolerance, 1e-6)))))
    strips: list[tuple[float, float]] = []
    for u in u_grid:
        interval = _v_range_at(polygon, u)
        if interval is None:
            return None
        strips.append(interval)
    triangles: list[tuple[Vec, Vec, Vec]] = []
    for index in range(len(u_grid) - 1):
        (a_low, a_high) = strips[index]
        (b_low, b_high) = strips[index + 1]
        u_a = u_grid[index]
        u_b = u_grid[index + 1]
        for j in range(steps_v):
            t0 = j / steps_v
            t1 = (j + 1) / steps_v
            a0 = (u_a, a_low + (a_high - a_low) * t0)
            b0 = (u_b, b_low + (b_high - b_low) * t0)
            b1 = (u_b, b_low + (b_high - b_low) * t1)
            a1 = (u_a, a_low + (a_high - a_low) * t1)
            triangles.append(
                (
                    surface.evaluate(*a0),
                    surface.evaluate(*b0),
                    surface.evaluate(*b1),
                )
            )
            triangles.append(
                (
                    surface.evaluate(*a0),
                    surface.evaluate(*b1),
                    surface.evaluate(*a1),
                )
            )
    return triangles


def _triangulate_face(
    surface: _Surface,
    loops: Sequence[Sequence[Vec]],
    *,
    angular_tolerance: float,
    max_depth: int,
) -> _FaceResult:
    result = _FaceResult(surface_kind=surface.kind)
    param_loops: list[list[tuple[float, float]]] = []
    for loop in loops:
        if len(loop) < 2:
            continue
        mapped = [surface.inverse(point) for point in loop]
        param_loops.append(_dedupe(_unwrap(mapped, surface)))
    param_loops = [loop for loop in param_loops if len(loop) >= 3]
    if not param_loops:
        result.ambiguous = "no boundary loop long enough to bound a region"
        return result
    areas = [abs(_signed_area(loop)) for loop in param_loops]
    outer_index = max(range(len(param_loops)), key=lambda index: areas[index])
    if areas[outer_index] <= 1e-14:
        # One closed loop on a doubly curved surface bounds two caps and the file
        # does not say which one the solid keeps: refuse rather than invent a side.
        result.ambiguous = (
            "boundary loop collapses to zero area in the surface's parameter space "
            "(closed-in-v surface trimmed to a single loop): the trimmed side is ambiguous"
        )
        return result
    outer = list(param_loops[outer_index])
    if _signed_area(outer) < 0:
        outer.reverse()
    holes: list[list[tuple[float, float]]] = []
    for index, loop in enumerate(param_loops):
        if index == outer_index:
            continue
        hole = list(loop)
        if _signed_area(hole) > 0:
            hole.reverse()
        holes.append(hole)
    if not holes:
        strip = _wrap_strip_triangles(surface, outer, angular_tolerance)
        if strip is not None:
            result.triangles.extend(strip)
            return result
    merged = _bridge_holes(outer, holes) if holes else outer
    base_triangles = _ear_clip(merged)
    if base_triangles is None:
        result.ambiguous = "ear clipping found no valid ear: the trim is not a simple polygon"
        return result
    # The triangulation has to cover the region exactly; a sum that misses is a
    # tessellator bug, and reporting it here beats reporting the volume it implies.
    region_area = abs(_signed_area(outer)) - sum(abs(_signed_area(hole)) for hole in holes)
    covered = sum(abs(_cross2(a, b, c)) for a, b, c in base_triangles) / 2.0
    if region_area > 0.0 and abs(covered - region_area) > 1e-6 * region_area:
        result.ambiguous = (
            f"triangulation covers {covered:.6g} of a {region_area:.6g} parameter region: "
            "overlapping ears, refusing to report a volume"
        )
        return result
    for triangle in base_triangles:
        pieces, exhausted = _refine(triangle, surface, angular_tolerance, max_depth)
        result.depth_exhausted += exhausted
        for piece in pieces:
            if abs(_cross2(piece[0], piece[1], piece[2])) <= 1e-18:
                continue
            result.triangles.append(tuple(surface.evaluate(u, v) for u, v in piece))  # type: ignore[arg-type]
    return result


# --------------------------------------------------------------------------- #
# boundary walk
# --------------------------------------------------------------------------- #
def _sample_edge(model: _Model, edge_id: str, orientation: bool, angular_step: float) -> list[Vec]:
    """Points along an edge in wire-traversal order.

    Two independent flags decide the order. ``ORIENTED_EDGE.orientation`` says whether
    the loop uses the edge in its own ``edge_start -> edge_end`` direction, and
    ``EDGE_CURVE.same_sense`` says whether the edge's direction follows increasing
    curve parameter. Crossing them gives both the traversal order and the parameter
    direction, which is what makes a seam edge pushed twice come back the way it went.
    """
    part = model.part
    entity_type = part.type(edge_id)
    if entity_type != "EDGE_CURVE":
        raise StepUnsupportedStep(edge_id, {entity_type: "edge type not implemented"})
    args = part.args(edge_id)
    start_point = model.vertex_point(part.refs(args[1])[0])
    end_point = model.vertex_point(part.refs(args[2])[0])
    curve = model.curve(part.refs(args[3])[0])
    same_sense = _flag(args[4])
    forward = orientation
    increasing = same_sense if orientation else not same_sense
    t_start = curve.inverse(start_point)
    t_end = curve.inverse(end_point)
    t_begin = t_start if forward else t_end
    t_finish = t_end if forward else t_start
    first_point = start_point if forward else end_point
    last_point = end_point if forward else start_point
    direction = 1.0 if increasing else -1.0
    if curve.periodic:
        span = (t_finish - t_begin) % (2.0 * math.pi) if increasing else (t_begin - t_finish) % (2.0 * math.pi)
        if span <= 1e-12:
            span = 2.0 * math.pi  # one vertex: the edge is the whole closed curve
        steps = min(720, max(2, int(math.ceil(span / max(angular_step, 1e-6)))))
        params = [t_begin + direction * span * index / steps for index in range(steps + 1)]
        points = [curve.evaluate(t) for t in params]
    else:
        span = abs(t_finish - t_begin)
        if span <= 1e-12 or isinstance(curve, LineCurve):
            return [first_point, last_point]
        steps = min(256, max(4, int(math.ceil(span / max(angular_step, 1e-6)))))
        params = [t_begin + direction * span * index / steps for index in range(steps + 1)]
        points = [curve.evaluate(t) for t in params]
    points[0] = first_point
    points[-1] = last_point
    return points


def _loop_points(model: _Model, loop_id: str, angular_step: float) -> list[Vec]:
    part = model.part
    entity_type = part.type(loop_id)
    if entity_type == "POLY_LOOP":
        return [model.point(ref) for ref in part.refs(part.args(loop_id)[1])]
    if entity_type != "EDGE_LOOP":
        raise StepUnsupportedStep(loop_id, {entity_type: "loop type not implemented"})
    points: list[Vec] = []
    for oriented_id in part.refs(part.args(loop_id)[1]):
        oriented_type = part.type(oriented_id)
        if oriented_type == "ORIENTED_EDGE":
            oriented_args = part.args(oriented_id)
            edge_id = part.refs(oriented_args[3])[0]
            orientation = _flag(oriented_args[4])
        elif oriented_type == "EDGE_CURVE":
            edge_id = oriented_id
            orientation = True
        else:
            raise StepUnsupportedStep(oriented_id, {oriented_type: "loop element not implemented"})
        segment = _sample_edge(model, edge_id, orientation, angular_step)
        if points and _norm(_sub(segment[0], points[-1])) <= 1e-9:
            points.extend(segment[1:])
        else:
            points.extend(segment)
    if len(points) > 2 and _norm(_sub(points[0], points[-1])) <= 1e-9:
        points.pop()
    return points


def _face_result(
    model: _Model,
    face_id: str,
    *,
    angular_step: float,
    angular_tolerance: float,
    max_depth: int,
) -> _FaceResult:
    part = model.part
    entity_type = part.type(face_id)
    if entity_type not in ("ADVANCED_FACE", "FACE_SURFACE"):
        raise StepUnsupportedStep(face_id, {entity_type: "face type not implemented"})
    args = part.args(face_id)
    surface = model.surface(part.refs(args[2])[0])
    same_sense = _flag(args[3])
    loops: list[list[Vec]] = []
    for bound_id in part.refs(args[1]):
        bound_args = part.args(bound_id)
        loop_id = part.refs(bound_args[1])[0]
        bound_orientation = _flag(bound_args[2])
        points = _loop_points(model, loop_id, angular_step)
        if not bound_orientation:
            points.reverse()
        loops.append(points)
    result = _triangulate_face(surface, loops, angular_tolerance=angular_tolerance, max_depth=max_depth)
    if not same_sense:
        result.triangles = [(a, c, b) for a, b, c in result.triangles]
    return result


# --------------------------------------------------------------------------- #
# assembly and measurement
# --------------------------------------------------------------------------- #
def _weld(
    triangles: Sequence[tuple[Vec, Vec, Vec]],
    patches: Sequence[str],
    tolerance: float,
) -> tuple[tuple[Vec, ...], tuple[tuple[int, int, int], ...], tuple[str, ...]]:
    """Merge vertices within ``tolerance``, checking the neighbouring cells too.

    Plain grid rounding is not enough: two faces reach the same rim point through
    different surface evaluations, the coordinates differ by ~1e-17 m, and if that
    point sits on a cell boundary the two copies land in different cells. The mesh
    then has a crack that looks like a modelling error. Searching the 27 surrounding
    cells and testing the actual distance removes that failure mode.
    """
    quantum = tolerance if tolerance > 0.0 else 1e-9
    cells: dict[tuple[int, int, int], list[int]] = defaultdict(list)
    vertices: list[Vec] = []
    faces: list[tuple[int, int, int]] = []
    patch_ids: list[str] = []
    for triangle, patch in zip(triangles, patches):
        indices: list[int] = []
        for point in triangle:
            key = (
                int(math.floor(point[0] / quantum)),
                int(math.floor(point[1] / quantum)),
                int(math.floor(point[2] / quantum)),
            )
            found = None
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    for dz in (-1, 0, 1):
                        for candidate in cells.get((key[0] + dx, key[1] + dy, key[2] + dz), ()):
                            other = vertices[candidate]
                            if (
                                abs(other[0] - point[0]) <= tolerance
                                and abs(other[1] - point[1]) <= tolerance
                                and abs(other[2] - point[2]) <= tolerance
                            ):
                                found = candidate
                                break
                        if found is not None:
                            break
                    if found is not None:
                        break
                if found is not None:
                    break
            if found is None:
                found = len(vertices)
                vertices.append(point)
                cells[key].append(found)
            indices.append(found)
        if len(set(indices)) < 3:
            continue
        faces.append((indices[0], indices[1], indices[2]))
        patch_ids.append(patch)
    return tuple(vertices), tuple(faces), tuple(patch_ids)


def _euler_characteristic(mesh: Mesh) -> int | None:
    if not mesh.faces:
        return None
    edges: set[tuple[int, int]] = set()
    for a, b, c in mesh.faces:
        for u, v in ((a, b), (b, c), (c, a)):
            edges.add((u, v) if u < v else (v, u))
    return len(mesh.vertices) - len(edges) + len(mesh.faces)


_SI_PREFIX_FACTORS = {
    "$": 1.0,
    ".DECI.": 1e1,
    ".CENTI.": 1e2,
    ".MILLI.": 1e3,
    ".MICRO.": 1e6,
    ".NANO.": 1e9,
}


def _declared_length_scale(part: Part21) -> tuple[float | None, str]:
    """Read the file's own LENGTH_UNIT as ``(file units per metre, provenance)``."""
    scales: dict[str, float] = {}
    for entity_id, (entity_type, args) in part.statements.items():
        raw = part.complex.get(entity_id, "")
        match = re.search(r"SI_UNIT\s*\(\s*(\$|\.[A-Z]+\.)\s*,\s*\.METRE\.\s*\)", raw or "")
        if match is not None:
            factor = _SI_PREFIX_FACTORS.get(match.group(1))
            if factor is not None:
                scales[entity_id] = factor
            continue
        if entity_type == "SI_UNIT" and len(args) >= 2 and str(args[1]).upper() == ".METRE.":
            factor = _SI_PREFIX_FACTORS.get(str(args[0]).upper() if args[0] else "$")
            if factor is not None:
                scales[entity_id] = factor
    # Conversion-based units (inch/foot) express one file unit as a multiple of a
    # unit that is already in ``scales``.
    for entity_id, (entity_type, args) in part.statements.items():
        raw = part.complex.get(entity_id, "")
        if "CONVERSION_BASED_UNIT" not in (raw or entity_type):
            continue
        match = re.search(r"CONVERSION_BASED_UNIT\s*\(\s*'([^']*)'\s*,\s*#(\d+)\s*\)", raw or "")
        if match is not None:
            measure_id = match.group(2)
        elif entity_type == "CONVERSION_BASED_UNIT":
            refs = part.refs(args[1])
            if not refs:
                continue
            measure_id = refs[0]
        else:
            continue
        measure_args = part.args(measure_id)
        factor = None
        for value in measure_args:
            if isinstance(value, Typed) and value.values and isinstance(value.values[0], float):
                factor = value.values[0]
                break
            if isinstance(value, float):
                factor = value
        if factor is None:
            continue
        base_refs = [ref for ref in part.refs(measure_args) if ref in scales]
        if not base_refs:
            continue
        # One file unit equals ``factor`` base units, so the number of file units per
        # metre is (base units per metre) / factor: an inch is 25.4 mm, i.e. 1000/25.4
        # = 39.37 inches per metre, not 25.4 * 1000.
        scales[entity_id] = scales[base_refs[0]] / factor
    used: str | None = None
    # GLOBAL_UNIT_ASSIGNED_CONTEXT is normally part of a complex instance, so it has
    # no simple entity type: read its unit list out of the raw complex arguments.
    contexts = list(part.of_type("GLOBAL_UNIT_ASSIGNED_CONTEXT"))
    for entity_id, raw in part.complex.items():
        if "GLOBAL_UNIT_ASSIGNED_CONTEXT" in raw:
            contexts.append(entity_id)
    for context_id in contexts:
        if part.type(context_id):
            candidates = part.refs(part.args(context_id)[1])
        else:
            match = re.search(r"GLOBAL_UNIT_ASSIGNED_CONTEXT\s*\(\s*\((.*?)\)\s*\)", part.complex.get(context_id, ""))
            candidates = [Ref(value) for value in re.findall(r"#(\d+)", match.group(1))] if match else []
        for ref in candidates:
            if ref in scales:
                used = ref
    if used is not None:
        return scales[used], f"file LENGTH_UNIT #{used}"
    if len(scales) == 1:
        only = next(iter(scales))
        return scales[only], f"file LENGTH_UNIT #{only}"
    if scales:
        return None, "several length units and no GLOBAL_UNIT_ASSIGNED_CONTEXT to pick one"
    return None, "no SI or conversion-based LENGTH_UNIT entity"# --------------------------------------------------------------------------- #
# public result
# --------------------------------------------------------------------------- #
@dataclass
class StepBrepImport:
    """What an analytic STEP read actually produced, including what it refused."""

    mesh: Mesh
    patch_ids: tuple[str, ...]
    source: str
    solid_count: int
    shell_count: int
    face_count: int
    volume_m3: float | None
    complete: bool
    length_scale: float
    length_scale_source: str
    surface_counts: dict[str, int]
    entity_counts: dict[str, int]
    unsupported: dict[str, str]
    ambiguous: dict[str, int]
    watertight: bool
    boundary_edges: int
    euler_characteristic: int | None
    metadata: dict[str, Any] = field(default_factory=dict)

    def as_dict(self) -> dict[str, Any]:
        return {
            "schema": STEP_BREP_SCHEMA,
            "source": self.source,
            "solid_count": self.solid_count,
            "shell_count": self.shell_count,
            "face_count": self.face_count,
            "triangles": len(self.mesh.faces),
            "vertices": len(self.mesh.vertices),
            "volume_m3": self.volume_m3,
            "complete": self.complete,
            "watertight": self.watertight,
            "boundary_edges": self.boundary_edges,
            "euler_characteristic": self.euler_characteristic,
            "length_scale": self.length_scale,
            "length_scale_source": self.length_scale_source,
            "surface_counts": dict(sorted(self.surface_counts.items())),
            "unsupported": dict(sorted(self.unsupported.items())),
            "ambiguous": dict(sorted(self.ambiguous.items())),
            "metadata": self.metadata,
        }

    def to_patch_geometry(self) -> Any:
        from .exchange import PatchGeometry

        return PatchGeometry(
            mesh=self.mesh,
            patch_ids=self.patch_ids,
            source=self.source,
            metadata={"schema": STEP_BREP_SCHEMA, **self.as_dict()},
        )


def _boundary_edge_count(mesh: Mesh) -> int:
    counts: dict[tuple[int, int], int] = defaultdict(int)
    for a, b, c in mesh.faces:
        for u, v in ((a, b), (b, c), (c, a)):
            counts[(u, v) if u < v else (v, u)] += 1
    return sum(1 for count in counts.values() if count == 1)


def read_step_brep(
    path: str | Path,
    *,
    length_scale: float | None = None,
    angular_step: float = 0.12,
    angular_tolerance: float = 0.10,
    max_depth: int = 8,
    weld_tolerance: float = 1e-7,
    allow_unsupported: bool = False,
) -> StepBrepImport:
    """Read an AP203/AP214 analytic B-rep file into a welded triangle mesh.

    ``length_scale`` is *file units per metre*; left as ``None`` it comes from the
    file's own ``LENGTH_UNIT`` and its provenance lands in
    ``length_scale_source``. ``angular_step`` bounds the sample spacing along a
    curved edge and ``angular_tolerance`` the angular span of a curved triangle,
    so both are accuracy knobs whose effect is measurable (see the docs: the
    volume converges quadratically in ``angular_tolerance``).

    ``volume_m3`` is reported only when every face was tessellated *and* the
    surface closed -- a partial or open mesh has no volume to report and says so
    with ``None`` instead of a number.
    """
    source = Path(path)
    text = source.read_text(encoding="utf-8", errors="replace")
    if "ISO-10303-21" not in text:
        raise StepParseError(f"{source} is not an ISO-10303-21 (STEP) file")
    if "END-ISO-10303-21" not in text:
        raise StepParseError(
            f"{source} is truncated: the exchange structure has no END-ISO-10303-21 "
            "terminator, so its data section cannot be trusted"
        )
    part = Part21(text)
    declared, unit_detail = _declared_length_scale(part)
    if length_scale is None:
        if declared is None:
            scale = 1000.0
            scale_source = f"assumed millimetre convention ({unit_detail})"
        else:
            scale = declared
            scale_source = f"declared {unit_detail}"
    else:
        scale = float(length_scale)
        scale_source = "caller override"
    if scale == 0.0:
        raise StepParseError("length_scale must be non-zero")
    model = _Model(part, scale)

    shells: list[tuple[str, bool]] = []
    solid_ids = part.of_type("MANIFOLD_SOLID_BREP", "BREP_WITH_VOIDS")
    for solid_id in solid_ids:
        args = part.args(solid_id)
        outer_refs = part.refs(args[1])
        if not outer_refs:
            continue
        shells.append((outer_refs[0], False))
        if part.type(solid_id) == "BREP_WITH_VOIDS" and len(args) > 2:
            for void_id in part.refs(args[2]):
                shells.append((void_id, True))
    if not shells:
        for shell_id in part.of_type("CLOSED_SHELL", "OPEN_SHELL"):
            shells.append((shell_id, False))
    if not shells:
        raise StepError(
            f"{source}: no MANIFOLD_SOLID_BREP / BREP_WITH_VOIDS / CLOSED_SHELL found; "
            "a wireframe or sheet-model file has no solid to measure"
        )

    triangles: list[tuple[Vec, Vec, Vec]] = []
    patches: list[str] = []
    unsupported: dict[str, str] = {}
    ambiguous: dict[str, int] = {}
    surface_counts: dict[str, int] = defaultdict(int)
    refused_faces = 0
    depth_exhausted = 0
    face_count = 0
    for shell_id, is_void in shells:
        shell_args = part.args(shell_id)
        for face_id in part.refs(shell_args[1]):
            face_count += 1
            try:
                outcome = _face_result(
                    model,
                    face_id,
                    angular_step=angular_step,
                    angular_tolerance=angular_tolerance,
                    max_depth=max_depth,
                )
            except StepUnsupportedStep as error:
                unsupported.update(error.unsupported)
                refused_faces += 1
                continue
            if outcome.ambiguous:
                ambiguous[outcome.ambiguous] = ambiguous.get(outcome.ambiguous, 0) + 1
                refused_faces += 1
                continue
            surface_counts[outcome.surface_kind] += 1
            depth_exhausted += outcome.depth_exhausted
            for triangle in outcome.triangles:
                if is_void:
                    triangle = (triangle[0], triangle[2], triangle[1])
                triangles.append(triangle)
                patches.append(face_id)

    if unsupported and not allow_unsupported:
        raise StepUnsupportedStep(str(source), unsupported)
    if ambiguous and not allow_unsupported:
        detail = "; ".join(f"{reason} (x{count})" for reason, count in ambiguous.items())
        raise StepTessellationError(f"{source}: ambiguous trims, no face invented: {detail}")

    vertices, faces, patch_ids = _weld(triangles, patches, weld_tolerance)
    mesh = Mesh(vertices=vertices, faces=faces)
    complete = refused_faces == 0
    watertight = bool(mesh.faces) and is_watertight(mesh)
    volume: float | None = None
    orientation_flipped = False
    if complete and watertight:
        volume = total_signed_volume(mesh)
        if volume < 0.0:
            mesh = Mesh(vertices=vertices, faces=tuple((a, c, b) for a, b, c in faces))
            volume = -volume
            orientation_flipped = True

    entity_counts: dict[str, int] = defaultdict(int)
    for entity_type, _ in part.statements.values():
        entity_counts[entity_type or "(COMPLEX)"] += 1
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    metadata = {
        "entities": len(part.statements),
        "header": part.header,
        "sha256": digest,
        "angular_step": angular_step,
        "angular_tolerance": angular_tolerance,
        "max_depth": max_depth,
        "weld_tolerance": weld_tolerance,
        "refused_faces": refused_faces,
        "depth_exhausted_triangles": depth_exhausted,
        "orientation_flipped": orientation_flipped,
        "volume_method": "divergence theorem over the tessellated boundary",
        "validation": "internal (analytic nominals + convergence); no third-party CAD kernel available",
    }
    return StepBrepImport(
        mesh=mesh,
        patch_ids=patch_ids,
        source=str(source),
        solid_count=len(solid_ids),
        shell_count=len(shells),
        face_count=face_count,
        volume_m3=volume,
        complete=complete,
        length_scale=scale,
        length_scale_source=scale_source,
        surface_counts=dict(surface_counts),
        entity_counts=dict(entity_counts),
        unsupported=dict(unsupported),
        ambiguous=dict(ambiguous),
        watertight=watertight,
        boundary_edges=_boundary_edge_count(mesh),
        euler_characteristic=_euler_characteristic(mesh),
        metadata=metadata,
    )


# --------------------------------------------------------------------------- #
# coverage census
# --------------------------------------------------------------------------- #
def census(path: str | Path) -> dict[str, Any]:
    """Entity census of a STEP file: what it contains and what this reader covers.

    Cheap (no tessellation) and honest: it reports every entity type in the file,
    which of them are supported surfaces/curves, and which ones this reader would
    refuse -- the machine-readable form of the "coverage list + unsupported types"
    acceptance criterion.
    """
    source = Path(path)
    text = source.read_text(encoding="utf-8", errors="replace")
    if "ISO-10303-21" not in text:
        raise StepParseError(f"{source} is not an ISO-10303-21 (STEP) file")
    part = Part21(text)
    counts: dict[str, int] = defaultdict(int)
    for entity_type, _ in part.statements.values():
        counts[entity_type or "(COMPLEX)"] += 1
    unsupported: dict[str, str] = {}
    for entity_type in counts:
        if entity_type in _SURFACE_REASONS:
            unsupported[entity_type] = _SURFACE_REASONS[entity_type]
        elif entity_type in _CURVE_REASONS:
            unsupported[entity_type] = _CURVE_REASONS[entity_type]
    declared, unit_detail = _declared_length_scale(part)
    surfaces = {name: counts[name] for name in SUPPORTED_SURFACES if counts.get(name)}
    curves = {name: counts[name] for name in SUPPORTED_CURVES if counts.get(name)}
    return {
        "source": str(source),
        "sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "bytes": source.stat().st_size,
        "header": part.header,
        "entities": dict(sorted(counts.items())),
        "instances": len(part.statements),
        "length_scale": declared,
        "length_scale_detail": unit_detail,
        "solids": {
            "MANIFOLD_SOLID_BREP": counts.get("MANIFOLD_SOLID_BREP", 0),
            "BREP_WITH_VOIDS": counts.get("BREP_WITH_VOIDS", 0),
            "SHELL_BASED_SURFACE_MODEL": counts.get("SHELL_BASED_SURFACE_MODEL", 0),
        },
        "shells": {name: counts.get(name, 0) for name in ("CLOSED_SHELL", "OPEN_SHELL")},
        "faces": {name: counts.get(name, 0) for name in ("ADVANCED_FACE", "FACE_SURFACE")},
        "supported_surfaces": surfaces,
        "supported_curves": curves,
        "unsupported": dict(sorted(unsupported.items())),
        "supported_surface_table": list(SUPPORTED_SURFACES),
        "supported_curve_table": list(SUPPORTED_CURVES),
    }


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def _write_obj(path: str | Path, mesh: Mesh) -> None:
    lines = ["# tensorcad.step_brep triangle mesh"]
    for x, y, z in mesh.vertices:
        lines.append(f"v {x!r} {y!r} {z!r}")
    for a, b, c in mesh.faces:
        lines.append(f"f {a + 1} {b + 1} {c + 1}")
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m tensorcad.step_brep",
        description="Read AP203/AP214 analytic B-rep STEP files; report coverage or mesh a solid.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    census_parser = subparsers.add_parser("census", help="entity census + supported/unsupported lists")
    census_parser.add_argument("path")
    mesh_parser = subparsers.add_parser("mesh", help="tessellate a solid and report its volume")
    mesh_parser.add_argument("path")
    mesh_parser.add_argument("--json", dest="json_path", default=None)
    mesh_parser.add_argument("--obj", dest="obj_path", default=None)
    mesh_parser.add_argument("--angular-step", type=float, default=0.12)
    mesh_parser.add_argument("--angular-tolerance", type=float, default=0.10)
    mesh_parser.add_argument("--allow-unsupported", action="store_true")
    args = parser.parse_args(argv)

    if args.command == "census":
        print(json.dumps(census(args.path), indent=2, sort_keys=True))
        return 0
    imported = read_step_brep(
        args.path,
        angular_step=args.angular_step,
        angular_tolerance=args.angular_tolerance,
        allow_unsupported=args.allow_unsupported,
    )
    payload = imported.as_dict()
    if args.json_path:
        Path(args.json_path).write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    if args.obj_path:
        _write_obj(args.obj_path, imported.mesh)
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())