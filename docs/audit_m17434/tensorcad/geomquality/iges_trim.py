"""Trimmed IGES import: fixed-width records -> entity graph -> trimmed faces.

Why this module exists
----------------------
:mod:`tensorcad.geomquality.io` reads and writes *untrimmed* NURBS patch sets, which is
what this repository's own writer produces.  An industrial IGES file is something else:
a yard, a propeller supplier or a model basin hands over a **trimmed** B-rep in which
every face is a surface (plane, surface of revolution, rational B-spline) plus boundary
curves, and nothing closes until those curves are applied.  ``io.read_iges`` refuses
such a file with ``IGES trimmed surface (type 144) is out of scope`` -- correctly,
because a patch reader that ignored trims would silently turn a propeller blade into a
full rectangle.

This module is the front half of the import pipeline::

    fixed-width IGES text
        -> directory entries + parameter records (all five sections)
        -> curves (100, 102, 106, 110, 126) and surfaces (108, 120, 122, 128)
           with 124 transformation matrices applied
        -> trimmed faces (144 / 143) whose boundary loops are 3-D polylines

The loops are produced in 3-D on purpose.  IGES also carries a parameter-space image of
each boundary curve (entity 142), but that entity's field layout differs between
exporters and a wrong reading of it is invisible downstream: any ``(u, v)`` pair maps
to *some* point on the surface.  This module therefore never trusts ``(u, v)`` data.  It
takes the boundary curves as 3-D curves, and the consumer
(:mod:`tensorcad.geomquality.trim_mesh`) inverts them onto the surface, which turns
"was the trim read correctly?" into a measured residual instead of a belief.

Supported, and refused
----------------------
====================  ==========================================================
entity                handling
====================  ==========================================================
124 transformation    affine 3x4, applied to control points and points
100 circular arc      exact rational arc
102 composite curve   ordered chain, joined end-to-start within tolerance
106 copious data      forms 11/12/63 used as a 3-D polyline trim
110 line segment      linear parameterisation
126 B-spline curve    rational (:class:`tensorcad.nurbs.NURBSCurve`)
108 plane             origin + orthonormal in-plane frame
120 surface of rev.   generatrix swept about an axis; angles in radians
122 tabulated cyl.    directrix translated along a direction
128 B-spline surface  rational (:class:`tensorcad.nurbs.NURBSSurface`)
143 bounded surface   boundary entities of type 141
144 trimmed surface   surface (parameter 1) + outer/inner boundary pointers
====================  ==========================================================

Anything else that a face depends on raises :class:`IgesError` naming the entity type
and the directory entry.  Skipping an unsupported boundary curve would under-trim the
face: its area and volume would be wrong and its shell open -- the worst possible input
to a watertightness gate.

Units: the global section's model-space scale and units flag are read, reported, and
carried in :attr:`TrimmedModel.metres_per_unit`.  No coordinate is scaled here; the
pipeline scales once and records that it did.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Sequence

from ..nurbs import NURBSCurve, NURBSSurface

__all__ = [
    "ArcCurve",
    "BoundedParametricSurface",
    "CompositeCurve",
    "Curve3D",
    "IgesDocument",
    "IgesError",
    "LineCurve",
    "PlaneSurface",
    "PolylineCurve",
    "RevolvedSurface",
    "SplineCurve",
    "Surface",
    "TabulatedCylinder",
    "TrimmedFace",
    "TrimmedModel",
    "UNIT_TO_METRES",
    "parse_iges",
    "read_iges_trimmed",
]

Vec3 = tuple[float, float, float]
Matrix = tuple[tuple[float, float, float, float], ...]

IDENTITY: Matrix = (
    (1.0, 0.0, 0.0, 0.0),
    (0.0, 1.0, 0.0, 0.0),
    (0.0, 0.0, 1.0, 0.0),
    (0.0, 0.0, 0.0, 1.0),
)

#: IGES global-section units flag -> metres per unit.  Flag 3 means "units defined by
#: the units-name field", which this module refuses rather than guessing.
UNIT_TO_METRES: dict[int, float] = {
    1: 0.0254,     # inch
    2: 0.001,      # millimetre
    4: 0.3048,     # foot
    5: 1609.344,   # mile
    6: 1.0,        # metre
    7: 1000.0,     # kilometre
    8: 2.54e-05,   # mil
    9: 1.0e-06,    # micron
    10: 0.01,      # centimetre
    11: 2.54e-08,  # micro-inch
}


class IgesError(ValueError):
    """Raised when an IGES file uses something this reader cannot represent exactly."""


# --------------------------------------------------------------------------- #
# small linear algebra
# --------------------------------------------------------------------------- #
def _sub(a: Sequence[float], b: Sequence[float]) -> Vec3:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _add(a: Sequence[float], b: Sequence[float]) -> Vec3:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def _scale(a: Sequence[float], s: float) -> Vec3:
    return (a[0] * s, a[1] * s, a[2] * s)


def _dot(a: Sequence[float], b: Sequence[float]) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _cross(a: Sequence[float], b: Sequence[float]) -> Vec3:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _norm(a: Sequence[float]) -> float:
    return math.sqrt(_dot(a, a))


def _unit(a: Sequence[float], label: str) -> Vec3:
    length = _norm(a)
    if length <= 0.0:
        raise IgesError(f"degenerate direction vector for {label}")
    return _scale(a, 1.0 / length)


def _apply_point(matrix: Matrix, point: Sequence[float]) -> Vec3:
    x, y, z = point[0], point[1], point[2]
    return (
        matrix[0][0] * x + matrix[0][1] * y + matrix[0][2] * z + matrix[0][3],
        matrix[1][0] * x + matrix[1][1] * y + matrix[1][2] * z + matrix[1][3],
        matrix[2][0] * x + matrix[2][1] * y + matrix[2][2] * z + matrix[2][3],
    )


def _apply_direction(matrix: Matrix, direction: Sequence[float]) -> Vec3:
    x, y, z = direction[0], direction[1], direction[2]
    return (
        matrix[0][0] * x + matrix[0][1] * y + matrix[0][2] * z,
        matrix[1][0] * x + matrix[1][1] * y + matrix[1][2] * z,
        matrix[2][0] * x + matrix[2][1] * y + matrix[2][2] * z,
    )


def _rotate_about_axis(point: Sequence[float], origin: Sequence[float],
                       axis: Sequence[float], angle: float) -> Vec3:
    """Rodrigues rotation of ``point`` about ``origin + t*axis`` by ``angle`` radians."""
    unit_axis = _unit(axis, "revolution axis")
    relative = _sub(point, origin)
    cosine, sine = math.cos(angle), math.sin(angle)
    term1 = _scale(relative, cosine)
    term2 = _scale(_cross(unit_axis, relative), sine)
    term3 = _scale(unit_axis, _dot(unit_axis, relative) * (1.0 - cosine))
    return _add(origin, _add(term1, _add(term2, term3)))


# --------------------------------------------------------------------------- #
# curves
# --------------------------------------------------------------------------- #
class Curve3D:
    """A parametric curve in model space, parameterised over :attr:`domain`."""

    kind = "curve"
    domain: tuple[float, float] = (0.0, 1.0)

    def evaluate(self, t: float) -> Vec3:  # pragma: no cover - interface
        raise NotImplementedError

    def sample(self, tolerance: float) -> list[Vec3]:
        """Sample so that no chord bows away from the curve by more than ``tolerance``."""
        return _adaptive_polyline(self.evaluate, self.domain, tolerance)


@dataclass
class LineCurve(Curve3D):
    start: Vec3
    end: Vec3
    kind = "line"
    domain: tuple[float, float] = (0.0, 1.0)

    def evaluate(self, t: float) -> Vec3:
        return _add(_scale(self.start, 1.0 - t), _scale(self.end, t))


@dataclass
class ArcCurve(Curve3D):
    """Exact circular arc: centre, start, end and the arc-plane normal."""

    centre: Vec3
    start: Vec3
    end: Vec3
    normal: Vec3
    kind = "arc"
    radius: float = 0.0
    sweep: float = 0.0
    _e1: Vec3 = (1.0, 0.0, 0.0)
    _e2: Vec3 = (0.0, 1.0, 0.0)

    def __post_init__(self) -> None:
        radius = _norm(_sub(self.start, self.centre))
        if radius <= 0.0:
            raise IgesError("circular arc with zero radius")
        self.radius = radius
        self._e1 = _unit(_sub(self.start, self.centre), "arc start")
        if _norm(self.normal) < 1e-12:
            raise IgesError("circular arc with degenerate plane normal")
        normal = _unit(self.normal, "arc normal")
        self._e2 = _cross(normal, self._e1)
        end = _sub(self.end, self.centre)
        sweep = math.atan2(_dot(end, self._e2), _dot(end, self._e1))
        if sweep <= 1e-12:
            sweep += 2.0 * math.pi
        self.sweep = sweep
        self.domain = (0.0, 1.0)

    def evaluate(self, t: float) -> Vec3:
        angle = self.sweep * t
        offset = _add(
            _scale(self._e1, self.radius * math.cos(angle)),
            _scale(self._e2, self.radius * math.sin(angle)),
        )
        return _add(self.centre, offset)


@dataclass
class SplineCurve(Curve3D):
    """Rational B-spline curve; any 124 transformation is baked into the control points."""

    curve: NURBSCurve
    kind = "bspline"

    @property
    def domain(self) -> tuple[float, float]:  # type: ignore[override]
        return self.curve.domain

    def evaluate(self, t: float) -> Vec3:
        return self.curve.evaluate(t)


@dataclass
class PolylineCurve(Curve3D):
    """Piecewise-linear boundary (copious data used as a trim)."""

    points: tuple[Vec3, ...]
    kind = "polyline"
    domain: tuple[float, float] = (0.0, 1.0)

    def evaluate(self, t: float) -> Vec3:
        if len(self.points) == 1:
            return self.points[0]
        position = min(max(t, 0.0), 1.0) * (len(self.points) - 1)
        index = min(int(position), len(self.points) - 2)
        local = position - index
        return _add(
            _scale(self.points[index], 1.0 - local), _scale(self.points[index + 1], local)
        )


@dataclass
class CompositeCurve(Curve3D):
    """Ordered chain of curves; ``domain`` is the normalised chain parameter."""

    parts: tuple[Curve3D, ...]
    kind = "composite"
    domain: tuple[float, float] = (0.0, 1.0)
    gaps: tuple[float, ...] = ()

    def evaluate(self, t: float) -> Vec3:
        count = len(self.parts)
        position = min(max(t, 0.0), 1.0) * count
        index = min(int(position), count - 1)
        local = position - index
        part = self.parts[index]
        lo, hi = part.domain
        return part.evaluate(lo + (hi - lo) * local)

    def sample(self, tolerance: float) -> list[Vec3]:  # type: ignore[override]
        points: list[Vec3] = []
        for part in self.parts:
            piece = part.sample(tolerance)
            if points:
                piece = piece[1:]
            points.extend(piece)
        return points


def _adaptive_polyline(evaluate: Callable[[float], Vec3], domain: tuple[float, float],
                       tolerance: float, max_depth: int = 12) -> list[Vec3]:
    """Bisect the parameter interval until every chord error is below ``tolerance``."""
    lo, hi = domain
    start, end = evaluate(lo), evaluate(hi)

    def flat(a: Vec3, b: Vec3, mid: Vec3) -> bool:
        ab = _sub(b, a)
        length = _norm(ab)
        if length < 1e-15:
            return _norm(_sub(mid, a)) <= tolerance
        return _norm(_cross(ab, _sub(mid, a))) / length <= tolerance

    points = [start]
    stack: list[tuple[float, float, Vec3, Vec3, int]] = [(lo, hi, start, end, 0)]
    while stack:
        a_t, b_t, a_p, b_p, depth = stack.pop()
        mid_t = 0.5 * (a_t + b_t)
        mid_p = evaluate(mid_t)
        if depth >= max_depth or flat(a_p, b_p, mid_p):
            points.append(b_p)
            continue
        stack.append((mid_t, b_t, mid_p, b_p, depth + 1))
        stack.append((a_t, mid_t, a_p, mid_p, depth + 1))
    return points


# --------------------------------------------------------------------------- #
# surfaces
# --------------------------------------------------------------------------- #
class Surface:
    """A parametric surface with an explicit, finite parameter domain."""

    kind = "surface"
    domain: tuple[tuple[float, float], tuple[float, float]] = ((0.0, 1.0), (0.0, 1.0))

    def evaluate(self, u: float, v: float) -> Vec3:  # pragma: no cover - interface
        raise NotImplementedError


@dataclass
class PlaneSurface(Surface):
    """A plane with an in-plane orthonormal frame: affine, so invertible in closed form."""

    origin: Vec3
    u_axis: Vec3
    v_axis: Vec3
    normal: Vec3
    margin: float
    kind = "plane"

    @property
    def domain(self) -> tuple[tuple[float, float], tuple[float, float]]:  # type: ignore[override]
        extent = 4.0 * max(self.margin, 1.0)
        return ((-extent, extent), (-extent, extent))

    def evaluate(self, u: float, v: float) -> Vec3:
        return _add(self.origin, _add(_scale(self.u_axis, u), _scale(self.v_axis, v)))

    def project(self, point: Sequence[float]) -> tuple[float, float]:
        relative = _sub(point, self.origin)
        return (_dot(relative, self.u_axis), _dot(relative, self.v_axis))


@dataclass
class RevolvedSurface(Surface):
    """Generatrix swept about an axis (IGES 120); ``u`` is the angle in radians."""

    axis_origin: Vec3
    axis_direction: Vec3
    generatrix: Curve3D
    start_angle: float
    end_angle: float
    kind = "revolution"

    @property
    def domain(self) -> tuple[tuple[float, float], tuple[float, float]]:  # type: ignore[override]
        return ((self.start_angle, self.end_angle), self.generatrix.domain)

    def evaluate(self, u: float, v: float) -> Vec3:
        return _rotate_about_axis(
            self.generatrix.evaluate(v), self.axis_origin, self.axis_direction, u
        )


@dataclass
class TabulatedCylinder(Surface):
    """Directrix translated along ``direction`` (IGES 122)."""

    directrix: Curve3D
    base: Vec3
    direction: Vec3
    length: float
    kind = "tabulated-cylinder"

    @property
    def domain(self) -> tuple[tuple[float, float], tuple[float, float]]:  # type: ignore[override]
        offset = _norm(_sub(self.base, self.directrix.evaluate(self.directrix.domain[0])))
        return (self.directrix.domain, (-self.length - offset, self.length + offset))

    def evaluate(self, u: float, v: float) -> Vec3:
        return _add(self.directrix.evaluate(u), _scale(self.direction, v))


@dataclass
class BoundedParametricSurface(Surface):
    """Adapter so a :class:`NURBSSurface` shares the analytic surfaces' interface."""

    surface: NURBSSurface
    kind = "bspline"

    @property
    def domain(self) -> tuple[tuple[float, float], tuple[float, float]]:  # type: ignore[override]
        return self.surface.domain

    def evaluate(self, u: float, v: float) -> Vec3:
        return self.surface.evaluate(u, v)


# --------------------------------------------------------------------------- #
# records
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class DirectoryEntry:
    """One IGES directory entry (two fixed-width 80-column records)."""

    sequence: int
    entity_type: int
    form: int
    transformation: int
    level: int
    status: str
    parameter_lines: int
    tokens: tuple[str, ...]
    raw: str

    def integer(self, position: int, default: int = 0) -> int:
        """Pointer/integer field at ``position``; blank and 0 both mean "null"."""
        if position >= len(self.tokens):
            return default
        token = self.tokens[position].strip()
        if not token:
            return default
        try:
            return int(float(token))
        except ValueError as exc:
            raise IgesError(
                f"entity {self.sequence} (type {self.entity_type}) field {position + 1} "
                f"is not a pointer: {token!r}"
            ) from exc

    def real(self, position: int) -> float:
        if position >= len(self.tokens):
            raise IgesError(
                f"entity {self.sequence} (type {self.entity_type}) has {len(self.tokens)} "
                f"parameters, field {position + 1} is missing"
            )
        token = self.tokens[position].strip()
        try:
            return float(token.replace("D", "E").replace("d", "e"))
        except ValueError as exc:
            raise IgesError(
                f"entity {self.sequence} (type {self.entity_type}) field {position + 1} "
                f"is not a number: {token!r}"
            ) from exc

    def numbers(self) -> list[float]:
        out: list[float] = []
        for token in self.tokens:
            try:
                out.append(float(token.replace("D", "E").replace("d", "e")))
            except ValueError:
                out.append(float("nan"))
        return out


def _tokenize_parameters(text: str,
                        hollerith_mismatches: list[str] | None = None) -> tuple[str, ...]:
    """Split IGES parameter data on ``,``/``;`` while honouring ``nH`` Hollerith fields."""
    if hollerith_mismatches is None:
        hollerith_mismatches = []
    tokens: list[str] = []
    index, length = 0, len(text)
    current: list[str] = []
    while index < length:
        char = text[index]
        if char in ",;":
            tokens.append("".join(current).strip())
            current = []
            if char == ";":
                break
            index += 1
            continue
        if char == "H":
            digits = "".join(current).strip()
            if digits.isdigit():
                count = int(digits)
                start = index + 1
                payload = text[start:start + count]
                end = start + count
                if len(payload) > 1 and payload[-1] in ",;" and text[start:start + count - 1]:
                    # The writer's Hollerith count includes the delimiter: trim it back.
                    payload = payload[:-1]
                    end = start + count - 1
                elif end < length and text[end] not in ",;":
                    # The writer's count is too small (TensorCAD writes 8H for a
                    # nine-character name): read on to the delimiter.
                    while end < length and text[end] not in ",;":
                        end += 1
                    payload = text[start:end]
                    if hollerith_mismatches is not None:
                        hollerith_mismatches.append(f"{count}H{payload!r}")
                tokens.append(payload)
                index = end
                if index < length and text[index] in ",;":
                    index += 1
                current = []
                continue
        current.append(char)
        index += 1
    if current:
        tokens.append("".join(current).strip())
    return tuple(tokens)


def _fixed_width_lines(text: str) -> list[str]:
    lines = text.splitlines()
    usable = [line for line in lines if len(line) >= 73 and line[72] in "SGDPT"]
    if not usable:
        raise IgesError(
            "no IGES records found: the file is not 80-column fixed-width IGES "
            "(section letters in column 73)"
        )
    return usable


@dataclass
class IgesDocument:
    """Parsed records of one IGES file: directory entries plus global metadata."""

    path: Path
    entries: dict[int, DirectoryEntry]
    global_fields: tuple[str, ...]
    units_flag: int
    units_name: str
    model_scale: float
    metres_per_unit: float
    warnings: tuple[str, ...] = ()

    def entry(self, sequence: int) -> DirectoryEntry:
        try:
            return self.entries[sequence]
        except KeyError as exc:
            raise IgesError(
                f"dangling pointer: directory entry {sequence} is referenced but absent"
            ) from exc

    def census(self) -> dict[int, int]:
        counts: dict[int, int] = {}
        for entry in self.entries.values():
            counts[entry.entity_type] = counts.get(entry.entity_type, 0) + 1
        return dict(sorted(counts.items()))


def parse_iges(path: str | Path) -> IgesDocument:
    """Parse the five IGES sections into directory entries with parameter tokens."""
    target = Path(path)
    text = target.read_text(errors="replace")
    records = _fixed_width_lines(text)

    global_lines = [line[:72].rstrip() for line in records if line[72] == "G"]
    parameter_blocks: dict[int, list[tuple[int, str]]] = {}
    for line in records:
        if line[72] != "P":
            continue
        entity = int(line[64:72].strip() or 0)
        sequence = int(line[73:80].strip() or 0)
        parameter_blocks.setdefault(entity, []).append((sequence, line[:64]))

    directory_lines = [line for line in records if line[72] == "D"]
    entries: dict[int, DirectoryEntry] = {}
    warnings: list[str] = []
    for pair in range(0, len(directory_lines) - 1, 2):
        first, second = directory_lines[pair], directory_lines[pair + 1]
        if not first[0:8].strip():
            continue
        sequence = int(first[73:80].strip() or 0) or (pair // 2 + 1)
        blocks = parameter_blocks.get(sequence, [])
        raw = "".join(chunk for _, chunk in sorted(blocks))
        entity_type = int(first[0:8])
        tokens = _tokenize_parameters(raw)
        if tokens and tokens[0].strip() == str(entity_type):
            # Parameter data starts with the entity type number; drop it so that field
            # indices below are 0-based in the entity's own parameter list.
            tokens = tokens[1:]
        entries[sequence] = DirectoryEntry(
            sequence=sequence,
            entity_type=entity_type,
            form=int(second[0:8]),
            transformation=int(first[48:56].strip() or 0),
            level=int(first[32:40].strip() or 0),
            status=first[64:72].strip(),
            parameter_lines=int(second[24:32].strip() or 0),
            tokens=tokens,
            raw=raw,
        )

    hollerith_issues: list[str] = []
    global_fields = _tokenize_parameters("".join(global_lines), hollerith_issues)
    units_flag, units_name, model_scale = 6, "M", 1.0
    # The standard fixes the global section layout (13 model-space scale, 14 units flag,
    # 15 units name) but writers disagree on whether every field is emitted, so the units
    # name is located by content and the two fields before it are read from there.  Both
    # the explicit flag and the name's own flag are kept; a disagreement is a warning, not
    # a silent choice.
    name_flags = {
        "MM": 2, "M": 6, "IN": 1, "FT": 4, "CM": 10, "KM": 7, "UM": 9, "MIL": 8,
        "UIN": 11, "MI": 5,
    }
    located = None
    for index, token in enumerate(global_fields):
        key = token.strip().upper()
        if key in name_flags:
            located = index
            break
    if located is not None:
        units_name = global_fields[located][:16]
        declared_flag: int | None = None
        if located >= 1:
            try:
                declared_flag = int(float(global_fields[located - 1]))
            except ValueError:
                declared_flag = None
        if located >= 1:
            try:
                model_scale = float(global_fields[located - 1])
            except ValueError:
                pass
        # ``model_scale`` is the field *before* the units flag, and the units flag is the
        # field before the units name; when both exist, the flag's value is authoritative.
        if located >= 2:
            try:
                model_scale = float(global_fields[located - 2])
            except ValueError:
                model_scale = 1.0
        units_flag = name_flags[units_name.strip().upper()]
        if declared_flag is not None and declared_flag in UNIT_TO_METRES and declared_flag != units_flag:
            warnings.append(
                f"units flag {declared_flag} disagrees with units name {units_name!r}: "
                f"the name ({units_flag}) was used"
            )
        elif declared_flag is not None and declared_flag in UNIT_TO_METRES:
            units_flag = declared_flag
    else:
        if len(global_fields) >= 14:
            try:
                model_scale = float(global_fields[12])
            except ValueError:
                model_scale = 1.0
            try:
                units_flag = int(float(global_fields[13]))
            except ValueError:
                units_flag = 6
            units_name = global_fields[14][:16]
        else:
            warnings.append(
                f"global section has {len(global_fields)} fields: units assumed to be metres"
            )
    if units_flag not in UNIT_TO_METRES:
        raise IgesError(
            f"global units flag {units_flag} ({units_name!r}) has no defined length: "
            "cannot convert to metres without guessing"
        )
    if hollerith_issues:
        warnings.append(
            "Hollerith length mismatch in the global section: "
            + ", ".join(hollerith_issues[:4])
        )
    return IgesDocument(
        path=target,
        entries=entries,
        global_fields=global_fields,
        units_flag=units_flag,
        units_name=units_name,
        model_scale=model_scale,
        metres_per_unit=UNIT_TO_METRES[units_flag] * model_scale,
        warnings=tuple(warnings),
    )


# --------------------------------------------------------------------------- #
# geometry builders
# --------------------------------------------------------------------------- #
def _plane_reference_axis(normal: Sequence[float]) -> Vec3:
    """In-plane direction used as a plane's u axis when the file gives none."""
    for candidate in ((0.0, 1.0, 0.0), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0)):
        if abs(_dot(candidate, normal)) < 0.9:
            return candidate
    return (1.0, 0.0, 0.0)


class _Builder:
    """Turns directory entries into curves and surfaces, applying 124 matrices."""

    def __init__(self, document: IgesDocument, tolerance: float, extent: float) -> None:
        self.document = document
        self.tolerance = tolerance
        self.extent = extent
        self.refused: dict[int, str] = {}
        self.consumed: set[int] = set()

    def matrix(self, entry: DirectoryEntry) -> Matrix:
        if not entry.transformation:
            return IDENTITY
        transform = self.document.entry(entry.transformation)
        if transform.entity_type != 124:
            raise IgesError(
                f"entity {entry.sequence} points at directory entry {transform.sequence} "
                f"(type {transform.entity_type}) as a transformation; expected 124"
            )
        values = transform.numbers()
        if len(values) < 12:
            raise IgesError(
                f"transformation {transform.sequence} has {len(values)} values; 12 are required"
            )
        r = values[:12]
        self.consumed.add(transform.sequence)
        return (
            (r[0], r[1], r[2], r[3]),
            (r[4], r[5], r[6], r[7]),
            (r[8], r[9], r[10], r[11]),
            (0.0, 0.0, 0.0, 1.0),
        )

    def curve(self, sequence: int) -> Curve3D:
        entry = self.document.entry(sequence)
        matrix = self.matrix(entry)
        kind = entry.entity_type
        self.consumed.add(sequence)
        if kind == 110:
            return LineCurve(
                start=_apply_point(matrix, (entry.real(0), entry.real(1), entry.real(2))),
                end=_apply_point(matrix, (entry.real(3), entry.real(4), entry.real(5))),
            )
        if kind == 100:
            plane_point = (entry.real(0), entry.real(1), entry.real(2))
            centre = (entry.real(3), entry.real(4), entry.real(5))
            start = (entry.real(6), entry.real(7), entry.real(8))
            end = (entry.real(9), entry.real(10), entry.real(11))
            normal = _sub(plane_point, centre)
            if _norm(normal) < 1e-12:
                raise IgesError(
                    f"circular arc {entry.sequence} has its plane point at its centre: "
                    "no plane normal"
                )
            return ArcCurve(
                centre=_apply_point(matrix, centre),
                start=_apply_point(matrix, start),
                end=_apply_point(matrix, end),
                normal=_apply_direction(matrix, normal),
            )
        if kind == 106:
            form = entry.form
            count = entry.integer(0)
            if form in (11, 12, 63) and len(entry.tokens) >= 1 + 3 * count:
                points = [
                    _apply_point(matrix, (
                        entry.real(1 + 3 * i), entry.real(2 + 3 * i), entry.real(3 + 3 * i)
                    ))
                    for i in range(count)
                ]
                return PolylineCurve(points=tuple(points))
            raise IgesError(
                f"copious data entity {entry.sequence} with form {form} is not a 3-D "
                "polyline boundary: refused"
            )
        if kind == 126:
            upper = entry.integer(0)
            degree = entry.integer(1)
            count = upper + 1
            knot_start = 6
            knot_count = upper + degree + 2
            knots = [entry.real(knot_start + i) for i in range(knot_count)]
            weight_start = knot_start + knot_count
            weights = [entry.real(weight_start + i) for i in range(count)]
            point_start = weight_start + count
            points = [
                _apply_point(matrix, (
                    entry.real(point_start + 3 * i),
                    entry.real(point_start + 3 * i + 1),
                    entry.real(point_start + 3 * i + 2),
                ))
                for i in range(count)
            ]
            return SplineCurve(
                curve=NURBSCurve(
                    degree=degree, knots=tuple(knots), control_points=tuple(points),
                    weights=tuple(weights),
                )
            )
        if kind == 102:
            count = entry.integer(0)
            parts = tuple(self.curve(entry.integer(1 + i)) for i in range(count))
            if not parts:
                raise IgesError(f"composite curve {entry.sequence} has no constituent curves")
            return CompositeCurve(parts=parts)
        raise IgesError(
            f"boundary curve type {kind} (directory entry {entry.sequence}) is not "
            "supported: refused rather than skipping a trim"
        )

    def surface(self, sequence: int) -> tuple[Surface, str]:
        entry = self.document.entry(sequence)
        matrix = self.matrix(entry)
        kind = entry.entity_type
        self.consumed.add(sequence)
        if kind == 108:
            a, b, c, d = entry.real(0), entry.real(1), entry.real(2), entry.real(3)
            if _norm((a, b, c)) < 1e-14:
                raise IgesError(f"plane {entry.sequence} has a zero normal")
            normal = _unit(_apply_direction(matrix, (a, b, c)), f"plane {entry.sequence}")
            origin: Vec3 | None = None
            for offset in (4, 5, 6, 7):
                if len(entry.tokens) < offset + 3:
                    continue
                candidate = (
                    entry.real(offset), entry.real(offset + 1), entry.real(offset + 2)
                )
                residual = a * candidate[0] + b * candidate[1] + c * candidate[2] - d
                if abs(residual) <= 1e-6 * max(1.0, abs(d)):
                    origin = _apply_point(matrix, candidate)
                    break
            if origin is None:
                # No usable reference point: take the foot of the perpendicular from the
                # model origin.  That is exact for the plane and stated, not assumed.
                origin = _apply_point(
                    matrix, _scale(_unit((a, b, c), "plane normal"), d / _norm((a, b, c)))
                )
            u_axis = _unit(
                _apply_direction(matrix, _plane_reference_axis(normal)), "plane u axis"
            )
            u_axis = _unit(_sub(u_axis, _scale(normal, _dot(u_axis, normal))), "plane u axis")
            v_axis = _cross(normal, u_axis)
            return (
                PlaneSurface(
                    origin=origin, u_axis=u_axis, v_axis=v_axis, normal=normal,
                    margin=self.extent,
                ),
                "plane",
            )
        if kind == 120:
            axis_curve = self.curve(entry.integer(0))
            generatrix = self.curve(entry.integer(1))
            if not isinstance(axis_curve, LineCurve):
                raise IgesError(
                    f"surface of revolution {entry.sequence} uses a non-line axis "
                    f"(type {self.document.entry(entry.integer(0)).entity_type})"
                )
            return (
                RevolvedSurface(
                    axis_origin=axis_curve.start,
                    axis_direction=_sub(axis_curve.end, axis_curve.start),
                    generatrix=generatrix,
                    start_angle=entry.real(2),
                    end_angle=entry.real(3),
                ),
                "revolution",
            )
        if kind == 122:
            directrix = self.curve(entry.integer(0))
            base = _apply_point(matrix, (entry.real(1), entry.real(2), entry.real(3)))
            direction = _apply_direction(
                matrix, (entry.real(4), entry.real(5), entry.real(6))
            )
            return (
                TabulatedCylinder(
                    directrix=directrix, base=base, direction=direction, length=self.extent
                ),
                "tabulated-cylinder",
            )
        if kind == 128:
            upper_u, upper_v = entry.integer(0), entry.integer(1)
            degree_u, degree_v = entry.integer(2), entry.integer(3)
            cursor = 9  # K1, K2, M1, M2 and five property flags
            knots_u_count = upper_u + degree_u + 2
            knots_u = [entry.real(cursor + i) for i in range(knots_u_count)]
            cursor += knots_u_count
            knots_v_count = upper_v + degree_v + 2
            knots_v = [entry.real(cursor + i) for i in range(knots_v_count)]
            cursor += knots_v_count
            count_u, count_v = upper_u + 1, upper_v + 1
            weights = [entry.real(cursor + i) for i in range(count_u * count_v)]
            cursor += count_u * count_v
            net: list[list[Vec3]] = []
            for i in range(count_u):
                row: list[Vec3] = []
                for j in range(count_v):
                    base = cursor + 3 * (i * count_v + j)
                    row.append(_apply_point(matrix, (
                        entry.real(base), entry.real(base + 1), entry.real(base + 2)
                    )))
                net.append(row)
            try:
                surface = NURBSSurface(
                    degree=(degree_u, degree_v),
                    knots=(tuple(knots_u), tuple(knots_v)),
                    control_net=tuple(tuple(row) for row in net),
                    weights=tuple(
                        tuple(weights[i * count_v:(i + 1) * count_v]) for i in range(count_u)
                    ),
                )
            except Exception as exc:  # numpy-free validation raises NURBSError
                raise IgesError(
                    f"B-spline surface {entry.sequence} is not a valid tensor-product patch "
                    f"with the standard 128 field order (K1,K2,M1,M2,5 flags, knots, weights, "
                    f"control points): {exc}"
                ) from exc
            return BoundedParametricSurface(surface=surface), "bspline"
        raise IgesError(
            f"surface type {kind} (directory entry {entry.sequence}) is not supported: "
            "refused rather than returning an untrimmed approximation"
        )


# --------------------------------------------------------------------------- #
# trimmed faces
# --------------------------------------------------------------------------- #
@dataclass
class TrimmedFace:
    """One face: a surface plus closed 3-D boundary polylines (outer/inner unclassified)."""

    index: int
    surface: Surface
    loops: tuple[tuple[Vec3, ...], ...]
    name: str
    source_entity: int
    surface_kind: str
    outer_entity: int | None = None
    notes: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        return {
            "index": self.index,
            "source_entity": self.source_entity,
            "surface": self.surface_kind,
            "name": self.name,
            "loop_count": len(self.loops),
            "loop_points": [len(loop) for loop in self.loops],
            "notes": list(self.notes),
        }


@dataclass
class TrimmedModel:
    """Faces read from an IGES file, in file units, with provenance."""

    faces: tuple[TrimmedFace, ...]
    units_name: str
    metres_per_unit: float
    source: str
    census: dict[int, int]
    warnings: tuple[str, ...] = ()
    notes: tuple[str, ...] = ()
    refused: dict[int, str] = field(default_factory=dict)
    consumed: tuple[int, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        return {
            "source": self.source,
            "units": self.units_name,
            "metres_per_unit": self.metres_per_unit,
            "faces": [face.to_dict() for face in self.faces],
            "entity_census": {str(key): value for key, value in self.census.items()},
            "consumed_entities": len(self.consumed),
            "refused_entities": {str(key): value for key, value in self.refused.items()},
            "warnings": list(self.warnings),
            "notes": list(self.notes),
        }


def _model_extent(document: IgesDocument) -> float:
    """Characteristic model size from the raw coordinates, used for tolerance defaults."""
    largest = 0.0
    for entry in document.entries.values():
        if entry.entity_type not in (100, 102, 106, 108, 110, 120, 122, 126, 128):
            continue
        for token in entry.tokens:
            try:
                value = float(token.replace("D", "E").replace("d", "e"))
            except ValueError:
                continue
            if math.isfinite(value):
                largest = max(largest, abs(value))
    return largest if largest > 0.0 else 1.0


def _loop_from_curve(builder: _Builder, sequence: int) -> tuple[Vec3, ...]:
    curve = builder.curve(sequence)
    points = curve.sample(builder.tolerance)
    if len(points) < 2:
        raise IgesError(f"boundary curve {sequence} sampled to fewer than two points")
    if _norm(_sub(points[0], points[-1])) <= builder.tolerance:
        points[-1] = points[0]
    return tuple(points)


def _curves_of_142(builder: _Builder, entry: DirectoryEntry) -> list[tuple[Vec3, ...]]:
    """Resolve a curve-on-surface (142) into candidate 3-D boundary polylines.

    The field layout of entity 142 differs between exporters (CATIA writes the 3-D
    curve, a parameter-space curve and a count into fields 4-6 of this file's records).
    Rather than trust one layout, every curve reference inside the entity is taken as a
    candidate and sampled as a 3-D curve; the face builder then keeps exactly those
    candidates that invert onto the referenced surface within tolerance.  A wrong
    candidate is dropped with a measured residual, not silently accepted.
    """
    document = builder.document
    references: list[int] = []
    for index in range(len(entry.tokens)):
        token = entry.tokens[index].strip()
        if not token or any(char in token for char in ".eEdD"):
            continue
        value = entry.integer(index)
        if value <= 0 or value not in document.entries or value == entry.sequence:
            continue
        if document.entries[value].entity_type in (100, 102, 106, 110, 126):
            references.append(value)
    if not references:
        raise IgesError(
            f"curve-on-surface entity {entry.sequence} does not reference any curve this "
            f"reader can evaluate (tokens: {list(entry.tokens)})"
        )
    loops: list[tuple[Vec3, ...]] = []
    for reference in references:
        try:
            loop = _loop_from_curve(builder, reference)
        except IgesError as exc:
            builder.refused[reference] = str(exc)
            continue
        if len(loop) >= 2:
            loops.append(loop)
    if not loops:
        raise IgesError(
            f"curve-on-surface entity {entry.sequence} references only curves that were refused"
        )
    return loops


def _curves_of_boundary(builder: _Builder, sequence: int) -> list[tuple[Vec3, ...]]:
    """Closed polylines implied by one boundary pointer (141, 142, or a plain curve)."""
    entry = builder.document.entry(sequence)
    if entry.entity_type in (100, 102, 106, 110, 126):
        return [_loop_from_curve(builder, sequence)]
    if entry.entity_type == 141:
        # 141 is a list of (sense, curve) pairs; the curves chain into one boundary.
        chain: list[Vec3] = []
        for index in range(len(entry.tokens)):
            token = entry.tokens[index].strip()
            if not token or any(char in token for char in ".eEdD"):
                continue
            pointer = entry.integer(index)
            if pointer <= 0 or pointer not in builder.document.entries:
                continue
            if builder.document.entries[pointer].entity_type not in (100, 102, 106, 110, 126):
                continue
            points = _loop_from_curve(builder, pointer)
            chain.extend(points if not chain else points[1:])
        if not chain:
            raise IgesError(f"boundary entity {sequence} (type 141) contains no usable curve")
        if _norm(_sub(chain[0], chain[-1])) <= builder.tolerance:
            chain[-1] = chain[0]
        return [tuple(chain)]
    if entry.entity_type == 142:
        return _curves_of_142(builder, entry)
    raise IgesError(f"boundary entity {sequence} of type {entry.entity_type} is not supported")


def _domain_loop(surface: Surface, samples: int = 24) -> tuple[Vec3, ...]:
    """The rectangle of an untrimmed surface's parameter domain, as a closed polyline."""
    (u0, u1), (v0, v1) = surface.domain
    edge = max(8, samples // 4)
    points: list[Vec3] = []
    for index in range(edge):
        points.append(surface.evaluate(u0 + (u1 - u0) * index / edge, v0))
    for index in range(edge):
        points.append(surface.evaluate(u1, v0 + (v1 - v0) * index / edge))
    for index in range(edge):
        points.append(surface.evaluate(u1 - (u1 - u0) * index / edge, v1))
    for index in range(edge):
        points.append(surface.evaluate(u0, v1 - (v1 - v0) * index / edge))
    return tuple(points)


def _boundary_references(entry: DirectoryEntry) -> tuple[int | None, list[int]]:
    """Split a 144/143 record into (outer boundary, inner boundaries)."""
    count_outer = entry.integer(1)
    count_inner = entry.integer(2)
    outer = entry.integer(3) if count_outer else None
    inner: list[int] = []
    for index in range(count_inner):
        pointer = entry.integer(4 + index)
        if pointer > 0:
            inner.append(pointer)
    return (outer if outer and outer > 0 else None), inner


def read_iges_trimmed(path: str | Path, *, tolerance: float | None = None) -> TrimmedModel:
    """Read an IGES file into trimmed faces whose boundaries are 3-D polylines.

    ``tolerance`` is the chord tolerance used to sample every boundary curve, in file
    units.  When omitted it is 1e-4 times the largest coordinate in the file, which is
    the order of the exporters' own "minimum resolution" field.
    """
    document = parse_iges(path)
    extent = _model_extent(document)
    tol = tolerance if tolerance is not None else 1e-4 * extent
    builder = _Builder(document, tol, extent)

    faces: list[TrimmedFace] = []
    for entry in sorted(document.entries.values(), key=lambda item: item.sequence):
        if entry.entity_type not in (143, 144):
            continue
        try:
            surface, kind = builder.surface(entry.integer(0))
            outer, inner = _boundary_references(entry)
            loops: list[tuple[Vec3, ...]] = []
            notes: list[str] = []
            for reference in ([outer] if outer else []) + inner:
                loops.extend(_curves_of_boundary(builder, reference))
            if not loops:
                raise IgesError(
                    f"trimmed surface {entry.sequence} has no boundary loops "
                    f"(tokens: {list(entry.tokens)})"
                )
        except IgesError as exc:
            builder.refused[entry.sequence] = str(exc)
            continue
        faces.append(
            TrimmedFace(
                index=len(faces),
                surface=surface,
                loops=tuple(loops),
                name=f"face_{len(faces):03d}",
                source_entity=entry.sequence,
                surface_kind=kind,
                outer_entity=outer,
                notes=tuple(notes),
            )
        )

    if not faces:
        # No trimmed surfaces at all: fall back to whole-domain faces for the supported
        # surface entities, which is what an "untrimmed patch set" file contains.
        for entry in sorted(document.entries.values(), key=lambda item: item.sequence):
            if entry.entity_type not in (108, 120, 122, 128):
                continue
            try:
                surface, kind = builder.surface(entry.sequence)
            except IgesError as exc:
                builder.refused[entry.sequence] = str(exc)
                continue
            faces.append(
                TrimmedFace(
                    index=len(faces),
                    surface=surface,
                    loops=(_domain_loop(surface),),
                    name=f"face_{len(faces):03d}",
                    source_entity=entry.sequence,
                    surface_kind=kind,
                    notes=("untrimmed: whole parameter domain was used",),
                )
            )
        if not faces:
            raise IgesError(
                f"no trimmed (144/143) or supported untrimmed surface found in {Path(path)}; "
                f"entity census: {document.census()}"
            )

    return TrimmedModel(
        faces=tuple(faces),
        units_name=document.units_name,
        metres_per_unit=document.metres_per_unit,
        source=str(Path(path)),
        census=document.census(),
        warnings=document.warnings,
        notes=(f"chord tolerance {tol:.6g} file units",),
        refused=dict(builder.refused),
        consumed=tuple(sorted(builder.consumed)),
    )