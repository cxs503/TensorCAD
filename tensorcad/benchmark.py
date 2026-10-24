"""One measurement contract for every geometry representation TensorCAD can produce.

The question this module answers is the one issue #19 asks: *given a geometry task, which
representation should an agent use, and what evidence supports that choice?* It answers it
by putting the representations through the **same** harness -- same rays, same reference,
same metric definitions, same gates -- and by refusing to invent a number for anything it
could not measure.

Four representation slots, taken from the issue
-----------------------------------------------

``analytic_field``
    The exact closed form. Its role is *reference*: the surface radius of the reference
    geometry is a formula, so the harness's own error is what this row measures (it comes
    out at machine precision, which is the point -- it proves the measurement chain).
``mesh_local_edit``
    A triangle mesh, edited the way a mesh editor would: move the vertices inside the
    edit region, leave everything else alone, no re-projection. This is the control that
    shows what discretisation alone costs.
``neural_field``
    A versioned ``.ngo`` bundle. This Python process cannot query it without the
    ``neural`` extra (torch), and the module does not pretend otherwise: with a bundle it
    imports the metrics that produced it, verifying the bundle's hashes first, and labels
    them ``imported`` with the producer named; without a bundle the row is
    ``not_supported`` with a command to run.
``ocp_step``
    The OpenCASCADE path. Sewing is the only thing this repository can do with a triangle
    soup, and it says so: the row reports a *faceted* round trip, never a feature B-rep.

Status is per metric, not per row
---------------------------------

``measured`` -- produced by this run, by code under test.
``imported`` -- read from an artifact this run verified but did not recompute.
``not_supported`` -- the capability needs a dependency or a kernel this environment
does not have; the reason names the command that would provide it.
``not_run`` -- available, deliberately skipped (a quick run, a size cap).
``failed`` -- attempted, and it failed; the failure is in the value's place.

Nothing is ever encoded as zero. A missing measurement is a missing measurement, and a
gate that reads a missing value fails rather than passing quietly.

The two cases
-------------

**Case A, sphere with a local dent.** Ground truth is a formula, the edit is local, and
the protected region is the far side of the sphere -- the exact task a neural local edit
is supposed to do well, and the task a naive mesh edit is expected to do cheaply. Rays are
a fixed Fibonacci set, so the sampling error is identical for every representation and
cannot favour one.

**Case B, the parametric hull.** A design loop drives three dimensions to a target
displacement and two shape ratios, then the result is exported and re-read. This is the
engineering case: it measures the loop, the mesh's validity, the analytic volume error and
the export round trip, and it measures the grid dependence by regenerating at twice the
resolution.

Only the standard library is used, so the benchmark runs in the same lane as the core.
Triangle ray casting is a Möller--Trumbore loop with a per-triangle angular prefilter:
`tier` in the rules, no acceleration structure, and no NumPy.
"""
from __future__ import annotations

import hashlib
import importlib
import json
import math
import platform
import statistics
import sys
import tempfile
import time
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field, replace
from functools import lru_cache
from pathlib import Path
from typing import Any

from . import mesh as meshlib
from .geometry import Mesh, MeshBuilder

__all__ = [
    "BENCHMARK_SCHEMA",
    "BenchmarkResult",
    "CaseResult",
    "Gate",
    "Metric",
    "RepresentationResult",
    "dent_radius",
    "fibonacci_directions",
    "hull_handoff_case",
    "run_benchmark",
    "sphere_dent_case",
]

#: Schema tag written into ``results.json``; bump on incompatible change.
BENCHMARK_SCHEMA = "tensorcad.benchmark/1"

#: Metric statuses. ``measured`` and ``imported`` carry a value; the rest never do.
MEASURED = "measured"
IMPORTED = "imported"
NOT_SUPPORTED = "not_supported"
NOT_RUN = "not_run"
FAILED = "failed"


@dataclass(frozen=True)
class Metric:
    """One number, with where it came from and what would make it pass.

    ``tolerance`` is a *relative* tolerance unless the unit is a rate or a count, which is
    stated per row by the case that defines the gate; ``passed`` is filled in by the gate
    evaluation so a reader never has to re-derive it.
    """

    value: float | None
    unit: str = ""
    status: str = MEASURED
    reference: str = ""
    tolerance: float | None = None
    note: str = ""
    passed: bool | None = None

    def __post_init__(self) -> None:
        if self.status not in (MEASURED, IMPORTED, NOT_SUPPORTED, NOT_RUN, FAILED):
            raise ValueError(f"unknown metric status {self.status!r}")
        if self.status in (MEASURED, IMPORTED):
            if self.value is None:
                raise ValueError(f"a {self.status} metric needs a value")
            if not math.isfinite(self.value):
                raise ValueError(
                    f"metric value {self.value!r} is not finite; use status={FAILED!r} and a note"
                )
        elif self.value is not None:
            raise ValueError(f"a {self.status} metric must not carry a value, got {self.value!r}")

    @property
    def available(self) -> bool:
        return self.status in (MEASURED, IMPORTED)

    def within(self, tolerance: float) -> Metric:
        """Fill in ``passed`` against a relative tolerance; missing values never pass."""
        if not self.available or self.value is None:
            return _with(self, passed=False)
        return _with(self, passed=abs(self.value) <= tolerance, tolerance=tolerance)

    def at_most(self, limit: float) -> Metric:
        """Fill in ``passed`` against an absolute limit in the metric's own unit."""
        if not self.available or self.value is None:
            return _with(self, passed=False)
        return _with(self, passed=self.value <= limit, tolerance=limit)

    def at_least(self, limit: float) -> Metric:
        if not self.available or self.value is None:
            return _with(self, passed=False)
        return _with(self, passed=self.value >= limit, tolerance=limit)

    def equals(self, expected: float, tolerance: float = 1e-12) -> Metric:
        """Exact-in-practice comparison, for counts and topology invariants."""
        if not self.available or self.value is None:
            return _with(self, passed=False)
        return _with(self, passed=abs(self.value - expected) <= tolerance, tolerance=tolerance)

    def to_dict(self) -> dict[str, Any]:
        return {
            "value": self.value,
            "unit": self.unit,
            "status": self.status,
            "reference": self.reference,
            "tolerance": self.tolerance,
            "note": self.note,
            "passed": self.passed,
        }


def _with(metric: Metric, **changes: Any) -> Metric:
    """A copy with the gate verdict (and its threshold) filled in."""
    updated = replace(metric, **changes)
    return updated


def _missing(status: str, reason: str, unit: str = "") -> Metric:
    return Metric(value=None, unit=unit, status=status, note=reason)


@dataclass(frozen=True)
class RepresentationResult:
    """What one representation achieved on one case."""

    identifier: str
    role: str
    status: str
    metrics: Mapping[str, Metric] = field(default_factory=dict)
    seconds: float = 0.0
    artifacts: Mapping[str, str] = field(default_factory=dict)
    notes: Sequence[str] = ()

    def gate_failures(self) -> tuple[str, ...]:
        return tuple(
            f"{self.identifier}.{name}" for name, metric in self.metrics.items() if metric.passed is False
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "identifier": self.identifier,
            "role": self.role,
            "status": self.status,
            "seconds": self.seconds,
            "artifacts": dict(self.artifacts),
            "notes": list(self.notes),
            "metrics": {name: metric.to_dict() for name, metric in self.metrics.items()},
        }


@dataclass(frozen=True)
class Gate:
    """A pass/fail judgement with the number that produced it."""

    name: str
    passed: bool
    detail: str = ""

    def to_dict(self) -> dict[str, Any]:
        return {"name": self.name, "passed": self.passed, "detail": self.detail}


@dataclass(frozen=True)
class CaseResult:
    """One benchmark case: a question, a reference, and the representations tried."""

    identifier: str
    question: str
    reference: str
    representations: tuple[RepresentationResult, ...]
    gates: tuple[Gate, ...] = ()
    notes: Sequence[str] = ()
    seconds: float = 0.0

    def by_id(self, identifier: str) -> RepresentationResult:
        for result in self.representations:
            if result.identifier == identifier:
                return result
        raise KeyError(f"{self.identifier} has no representation {identifier!r}")

    @property
    def failed_gates(self) -> tuple[Gate, ...]:
        return tuple(gate for gate in self.gates if not gate.passed)

    def to_dict(self) -> dict[str, Any]:
        return {
            "identifier": self.identifier,
            "question": self.question,
            "reference": self.reference,
            "seconds": self.seconds,
            "notes": list(self.notes),
            "gates": [gate.to_dict() for gate in self.gates],
            "representations": [result.to_dict() for result in self.representations],
        }


def _environment() -> dict[str, Any]:
    """Versions and platform, so a number can be reproduced or explained away."""
    from . import __version__

    def version(module: str) -> str | None:
        try:
            loaded = __import__(module)
        except Exception:  # pragma: no cover - depends on the install
            return None
        return getattr(loaded, "__version__", "installed")

    return {
        "tensorcad": __version__,
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "machine": platform.machine(),
        "processor": platform.processor() or "unknown",
        "numpy": version("numpy"),
        "torch": version("torch"),
        "safetensors": version("safetensors"),
        "ocp": version("OCP"),
        "cpu_count": __import__("os").cpu_count(),
    }


@dataclass(frozen=True)
class BenchmarkResult:
    """Every case, every gate, and the environment they were measured in."""

    cases: tuple[CaseResult, ...]
    environment: Mapping[str, Any]
    seconds: float
    command: str = ""

    @property
    def passed(self) -> bool:
        return all(not case.failed_gates for case in self.cases)

    @property
    def failures(self) -> tuple[str, ...]:
        return tuple(f"{case.identifier}: {gate.name}" for case in self.cases for gate in case.failed_gates)

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": BENCHMARK_SCHEMA,
            "passed": self.passed,
            "failures": list(self.failures),
            "seconds": self.seconds,
            "command": self.command,
            "environment": dict(self.environment),
            "cases": [case.to_dict() for case in self.cases],
        }

    def write_json(self, path: str | Path) -> Path:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(self.to_dict(), indent=2) + "\n", encoding="utf-8")
        return target

    def write_markdown(self, path: str | Path) -> Path:
        """A generated report: same numbers as the JSON, no hand editing."""
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        lines: list[str] = [
            "# Geometry representation benchmark",
            "",
            "Generated from `results.json` by `tensorcad.benchmark`; every number below was",
            "produced by the run that wrote this file.",
            "",
            f"- command: `{self.command}`",
            f"- tensorcad {self.environment['tensorcad']}, python {self.environment['python']}, "
            f"numpy {self.environment['numpy']}, torch {self.environment['torch']}, "
            f"ocp {self.environment['ocp']}",
            f"- wall clock: {self.seconds:.2f} s",
            f"- **verdict: {'pass' if self.passed else 'FAIL'}**"
            + ("" if self.passed else f" — {len(self.failures)} gate(s) failed"),
            "",
        ]
        for case in self.cases:
            lines.extend(
                [
                    f"## {case.identifier}",
                    "",
                    f"*{case.question}*",
                    "",
                    f"Reference: {case.reference}",
                    "",
                ]
            )
            if case.notes:
                lines.extend([f"- {note}" for note in case.notes] + [""])
            for result in case.representations:
                lines.extend(
                    [
                        f"### `{result.identifier}` — {result.role} "
                        f"({result.status}, {result.seconds:.2f} s)",
                        "",
                        "| metric | value | unit | status | gate |",
                        "|---|---|---|---|---|",
                    ]
                )
                for name, metric in result.metrics.items():
                    value = "—" if metric.value is None else _format(metric.value)
                    gate = "—" if metric.passed is None else ("pass" if metric.passed else "**fail**")
                    lines.append(f"| `{name}` | {value} | {metric.unit} | {metric.status} | {gate} |")
                lines.append("")
                for note in result.notes:
                    lines.append(f"> {note}")
                if result.notes:
                    lines.append("")
            lines.extend(["#### gates", "", "| gate | result | detail |", "|---|---|---|"])
            for case_gate in case.gates:
                verdict = "pass" if case_gate.passed else "**fail**"
                lines.append(f"| `{case_gate.name}` | {verdict} | {case_gate.detail} |")
            lines.append("")
        target.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return target


def _format(value: float) -> str:
    if value == 0.0:
        return "0"
    if abs(value) < 1e-3 or abs(value) >= 1e5:
        return f"{value:.3e}"
    return f"{value:.6g}"


# --------------------------------------------------------------------------- #
# the measuring instruments
# --------------------------------------------------------------------------- #


def fibonacci_directions(count: int) -> tuple[tuple[float, float, float], ...]:
    """A deterministic, near-uniform set of unit directions (no RNG, no seed to lose).

    Every representation is probed along exactly these rays, so the sampling error is
    common to all of them and cannot favour one. The golden-angle spiral has a maximal
    minimum spacing for its size, which matters here: the dent subtends a small solid
    angle, and a poor sampling would under-resolve it for every row at once.
    """
    if count < 1:
        raise ValueError(f"count must be at least 1, got {count!r}")
    golden = math.pi * (3.0 - math.sqrt(5.0))
    directions: list[tuple[float, float, float]] = []
    for index in range(count):
        z = 1.0 - 2.0 * (index + 0.5) / count
        radius = math.sqrt(max(0.0, 1.0 - z * z))
        theta = golden * index
        directions.append((radius * math.cos(theta), radius * math.sin(theta), z))
    return tuple(directions)


def stratified_directions(count: int, cone_angle: float) -> tuple[tuple[float, float, float], ...]:
    """Half the rays uniform on the sphere, half spread inside ``cone_angle`` of +x.

    A uniform set is the fair way to probe the whole surface, but the edit subtends a small
    solid angle: at 512 uniform rays only four land inside the dent, which is too few to
    measure it. Stratifying fixes the resolution of the region under test without giving any
    representation an advantage, because every row is probed along the same rays.
    """
    if count < 4:
        raise ValueError(f"count must be at least 4, got {count!r}")
    if not 0.0 < cone_angle < math.pi:
        raise ValueError(f"cone_angle must be in (0, pi), got {cone_angle!r}")
    inside = count // 3
    outside = count - inside
    directions = list(fibonacci_directions(outside))
    # Equal-area spiral inside the spherical cap, so the cone is sampled without clumping
    # at its apex.
    golden = math.pi * (3.0 - math.sqrt(5.0))
    cos_edge = math.cos(cone_angle)
    for index in range(inside):
        cosine = 1.0 - (1.0 - cos_edge) * (index + 0.5) / inside
        sine = math.sqrt(max(0.0, 1.0 - cosine * cosine))
        azimuth = golden * index
        directions.append((cosine, sine * math.cos(azimuth), sine * math.sin(azimuth)))
    return tuple(directions)


def dent_radius(
    direction: Sequence[float], radius: float, depth: float, sigma: float
) -> float:
    """Radius of the reference surface along ``direction``: a smooth inward dent on +x.

    ``r(theta) = R - depth * exp(-0.5 * (theta / sigma)^2)``, with ``theta`` the angle from
    the +x axis and ``r`` the distance from the origin. This is the same reference used by
    ``examples/neural_geometry/edit_sphere_dent.py``, so a neural result and a mesh result
    are compared against the same target rather than each against its own.
    """
    length = math.sqrt(sum(component * component for component in direction))
    if length == 0.0:
        raise ValueError("direction must not be the zero vector")
    cosine = max(-1.0, min(1.0, direction[0] / length))
    angle = math.acos(cosine)
    return radius - depth * math.exp(-0.5 * (angle / sigma) ** 2)


def _dent_derivative(direction: Sequence[float], radius: float, depth: float, sigma: float) -> float:
    length = math.sqrt(sum(component * component for component in direction))
    cosine = max(-1.0, min(1.0, direction[0] / length))
    angle = math.acos(cosine)
    return depth * math.exp(-0.5 * (angle / sigma) ** 2) * angle / (sigma * sigma)


def _exact_volume(radius: float, depth: float, sigma: float, nodes: int = 4001) -> float:
    """``(1/3) * 2*pi * integral r(theta)^3 sin(theta) dtheta`` by Simpson's rule.

    The reference integral, not a mesh measurement: this is what every mesh volume in
    case A is compared against, so it is computed from the formula at a resolution where
    its own error is far below anything being measured (verified by refining in the tests).
    """
    total = 0.0
    for index in range(nodes):
        theta = math.pi * index / (nodes - 1)
        weight = 1.0 if index in (0, nodes - 1) else (4.0 if index % 2 else 2.0)
        r = radius - depth * math.exp(-0.5 * (theta / sigma) ** 2)
        total += weight * r**3 * math.sin(theta)
    return (2.0 * math.pi / 3.0) * total * (math.pi / (nodes - 1)) / 3.0


def _exact_area(radius: float, depth: float, sigma: float, nodes: int = 4001) -> float:
    """``2*pi * integral r sin(theta) sqrt(r^2 + (dr/dtheta)^2) dtheta``, via Simpson."""
    total = 0.0
    for index in range(nodes):
        theta = math.pi * index / (nodes - 1)
        weight = 1.0 if index in (0, nodes - 1) else (4.0 if index % 2 else 2.0)
        r = radius - depth * math.exp(-0.5 * (theta / sigma) ** 2)
        dr = depth * math.exp(-0.5 * (theta / sigma) ** 2) * theta / (sigma * sigma)
        total += weight * r * math.sqrt(r * r + dr * dr) * math.sin(theta)
    return 2.0 * math.pi * total * (math.pi / (nodes - 1)) / 3.0


class _RayCaster:
    """Nearest positive ray--triangle intersection from the origin, with a cheap prefilter.

    The surface here is star-shaped about the origin (a sphere with a shallow inward dent),
    so a ray crosses it exactly once and the nearest hit is the surface. The prefilter is
    per-triangle: the angle between the ray and the triangle's centroid must be within the
    triangle's own angular radius. That skips the far side of the sphere for free, which is
    most of the work at these ray counts; it is not an acceleration structure, and the
    class says so rather than implying one.
    """

    def __init__(self, mesh: Mesh) -> None:
        self.mesh = mesh
        triangles: list[tuple[tuple[float, float, float], ...]] = []
        centres: list[tuple[float, float, float]] = []
        cosines: list[float] = []
        for face in mesh.faces:
            a, b, c = (mesh.vertices[index] for index in face)
            triangles.append((a, b, c))
            centre = ((a[0] + b[0] + c[0]) / 3.0, (a[1] + b[1] + c[1]) / 3.0, (a[2] + b[2] + c[2]) / 3.0)
            distance = math.sqrt(sum(component * component for component in centre))
            if distance == 0.0:  # pragma: no cover - a degenerate face at the origin
                centres.append((0.0, 0.0, 1.0))
                cosines.append(-1.0)
                continue
            centres.append((centre[0] / distance, centre[1] / distance, centre[2] / distance))
            spread = max(
                math.dist(centre, vertex) for vertex in (a, b, c)
            ) / distance
            cosines.append(math.cos(min(math.pi, math.asin(min(1.0, spread)) + 0.05)))
        self._triangles = triangles
        self._centres = centres
        self._cosines = cosines

    def radius(self, direction: Sequence[float]) -> float | None:
        """Distance to the first crossing along ``direction``, or ``None`` if it misses."""
        nearest: float | None = None
        dx, dy, dz = direction
        for index, (a, b, c) in enumerate(self._triangles):
            centre = self._centres[index]
            if dx * centre[0] + dy * centre[1] + dz * centre[2] < self._cosines[index]:
                continue
            hit = _triangle_hit(a, b, c, dx, dy, dz)
            if hit is not None and (nearest is None or hit < nearest):
                nearest = hit
        return nearest


def _triangle_hit(
    a: Sequence[float],
    b: Sequence[float],
    c: Sequence[float],
    dx: float,
    dy: float,
    dz: float,
) -> float | None:
    """Möller--Trumbore for a ray from the origin; returns the parameter ``t``."""
    edge1 = (b[0] - a[0], b[1] - a[1], b[2] - a[2])
    edge2 = (c[0] - a[0], c[1] - a[1], c[2] - a[2])
    px = dy * edge2[2] - dz * edge2[1]
    py = dz * edge2[0] - dx * edge2[2]
    pz = dx * edge2[1] - dy * edge2[0]
    determinant = edge1[0] * px + edge1[1] * py + edge1[2] * pz
    if abs(determinant) < 1e-18:
        return None
    inverse = 1.0 / determinant
    tx = -a[0]
    ty = -a[1]
    tz = -a[2]
    u = (tx * px + ty * py + tz * pz) * inverse
    if u < -1e-12 or u > 1.0 + 1e-12:
        return None
    qx = ty * edge1[2] - tz * edge1[1]
    qy = tz * edge1[0] - tx * edge1[2]
    qz = tx * edge1[1] - ty * edge1[0]
    v = (dx * qx + dy * qy + dz * qz) * inverse
    if v < -1e-12 or u + v > 1.0 + 1e-12:
        return None
    t = (edge2[0] * qx + edge2[1] * qy + edge2[2] * qz) * inverse
    return t if t > 0.0 else None


def _lat_long_surface(
    radius_of: Callable[[tuple[float, float, float]], float],
    *,
    rings: int,
    segments: int,
) -> Mesh:
    """Triangulate a star-shaped surface sampled on a sphere of directions.

    ``radius_of`` is called with a unit direction and returns the surface radius there.
    Rings run pole to pole about +x (the dent axis), segments around it. Both poles are
    single vertices, so every edge is shared by exactly two triangles and the result is
    closed by construction -- no welding pass, no tolerance.
    """
    if rings < 2 or segments < 3:
        raise ValueError(f"need rings >= 2 and segments >= 3, got {rings} x {segments}")
    builder = MeshBuilder(units="m", coordinate_frame="object_local")

    def point(theta_index: int, segment: int) -> int:
        theta = math.pi * theta_index / rings
        phi = 2.0 * math.pi * segment / segments
        direction = (
            math.cos(theta),
            math.sin(theta) * math.cos(phi),
            math.sin(theta) * math.sin(phi),
        )
        r = radius_of(direction)
        return builder.add_vertex((r * direction[0], r * direction[1], r * direction[2]))

    top = point(0, 0)
    ring_index: list[list[int]] = []
    for theta_index in range(1, rings):
        ring_index.append([point(theta_index, segment) for segment in range(segments)])
    bottom = point(rings, 0)

    for segment in range(segments):
        following = (segment + 1) % segments
        builder.add_triangle(top, ring_index[0][segment], ring_index[0][following])
    for theta_index in range(len(ring_index) - 1):
        upper = ring_index[theta_index]
        lower = ring_index[theta_index + 1]
        for segment in range(segments):
            following = (segment + 1) % segments
            builder.add_quad(upper[segment], lower[segment], lower[following], upper[following])
    last = ring_index[-1]
    for segment in range(segments):
        following = (segment + 1) % segments
        builder.add_triangle(bottom, last[following], last[segment])

    mesh = builder.build()
    if meshlib.total_signed_volume(mesh) < 0.0:  # pragma: no cover - winding guard
        mesh = meshlib.orient_outward(mesh).mesh
    return mesh


def _digest(mesh: Mesh) -> str:
    import struct

    digest = hashlib.sha256()
    digest.update(struct.pack("<II", len(mesh.vertices), len(mesh.faces)))
    for vertex in mesh.vertices:
        digest.update(struct.pack("<3d", *vertex))
    for face in mesh.faces:
        digest.update(struct.pack("<3I", *face))
    return digest.hexdigest()


# --------------------------------------------------------------------------- #
# Case A: a local edit with a closed-form reference
# --------------------------------------------------------------------------- #

DENT_RADIUS_M = 0.5
DENT_DEPTH_M = 0.06
DENT_SIGMA_RAD = 0.22
#: Rays inside this cone are "the edit"; rays outside ``PROTECTED_FACTOR * sigma`` are
#: "the protected region". The gap between them is deliberately unjudged: the reference
#: dent has a Gaussian tail that never reaches zero, so a hard boundary at 1 sigma would
#: be a boundary the reference itself does not have.
DENT_EDGE_FACTOR = 1.0
PROTECTED_FACTOR = 2.5


def _angle_from_axis(direction: Sequence[float]) -> float:
    return math.acos(max(-1.0, min(1.0, direction[0])))


def sphere_dent_case(
    *,
    rays: int = 512,
    rings: int = 80,
    segments: int = 160,
    neural_bundle: str | Path | None = None,
    ocp_artifact: str | Path | None = None,
) -> CaseResult:
    """Case A: apply a known local dent, then measure who reproduced it and who did not.

    The four rows share the ray set, the reference formula and the metric definitions. The
    analytic row is the reference; the mesh row performs the edit the only way a mesh can
    (move the vertices inside the edit region); the neural row imports a verified bundle's
    own validation numbers when one is supplied, and says ``not_supported`` when none is;
    the OCP row measures a faceted STEP round trip only when the kernel is installed.
    """
    started = time.perf_counter()
    directions = stratified_directions(rays, DENT_EDGE_FACTOR * DENT_SIGMA_RAD)
    angles = [_angle_from_axis(direction) for direction in directions]
    dent_mask = [angle <= DENT_EDGE_FACTOR * DENT_SIGMA_RAD for angle in angles]
    # The edit has a hard support edge at PROTECTED_FACTOR * sigma, but a triangulated
    # surface does not: a triangle that straddles the edge carries the edit up to one cell
    # further out. A ray inside that band is not evidence of a leaky edit, it is evidence
    # that the edit boundary is a staircase. So the protected region starts one and a half
    # cells beyond the support edge (no reachable triangle then has a moved vertex), and the
    # band itself is measured and reported separately.
    cell_angle = max(math.pi / rings, 2.0 * math.pi / segments)
    protected_edge = PROTECTED_FACTOR * DENT_SIGMA_RAD
    protected_start = protected_edge + 1.5 * cell_angle
    protected_mask = [angle >= protected_start for angle in angles]
    transition_mask = [protected_edge <= angle < protected_start for angle in angles]
    if not any(dent_mask) or not any(protected_mask):  # pragma: no cover - needs a tiny ray count
        raise ValueError("the ray set does not resolve both the edit and the protected region")

    reference_radii = [
        dent_radius(direction, DENT_RADIUS_M, DENT_DEPTH_M, DENT_SIGMA_RAD)
        for direction in directions
    ]
    exact_volume = _exact_volume(DENT_RADIUS_M, DENT_DEPTH_M, DENT_SIGMA_RAD)
    exact_area = _exact_area(DENT_RADIUS_M, DENT_DEPTH_M, DENT_SIGMA_RAD)
    representations: list[RepresentationResult] = []

    # -- analytic field: the reference row ---------------------------------- #
    analytic_started = time.perf_counter()
    analytic_mesh = _lat_long_surface(
        lambda direction: dent_radius(direction, DENT_RADIUS_M, DENT_DEPTH_M, DENT_SIGMA_RAD),
        rings=rings,
        segments=segments,
    )
    analytic_seconds = time.perf_counter() - analytic_started
    analytic_volume = meshlib.total_signed_volume(analytic_mesh)
    analytic_area = meshlib.surface_area(analytic_mesh)
    analytic_metrics = {
        "reference_radius_error_max_m": Metric(
            0.0,
            "m",
            reference="the same formula the reference uses",
            note="ground truth by construction: this row proves the harness is consistent, "
            "it does not test an approximation",
        ),
        "volume_m3": Metric(exact_volume, "m^3", reference="quadrature of the formula", tolerance=1e-6),
        "surface_area_m2": Metric(exact_area, "m^2", reference="quadrature of the formula", tolerance=1e-6),
        "triangulated_volume_error": Metric(
            abs(analytic_volume / exact_volume - 1.0),
            "relative",
            reference="exact volume",
            note="an exact surface still has to be triangulated to be exported; this is the "
            f"cost at {rings}x{segments} rings/segments, the ideal parametrisation, not marching tets",
        ),
        "triangulated_area_error": Metric(
            abs(analytic_area / exact_area - 1.0), "relative", reference="exact area"
        ),
        "watertight": Metric(
            1.0 if meshlib.is_watertight(analytic_mesh) else 0.0, "bool", reference="every edge shared twice"
        ),
    }
    analytic_metrics["reference_radius_error_max_m"] = analytic_metrics[
        "reference_radius_error_max_m"
    ].at_most(1e-12)
    analytic_metrics["watertight"] = analytic_metrics["watertight"].at_least(1.0)
    representations.append(
        RepresentationResult(
            identifier="analytic_field",
            role="reference: exact closed form",
            status=MEASURED,
            metrics=analytic_metrics,
            seconds=analytic_seconds,
            artifacts={"triangulated_surface": _digest(analytic_mesh)},
            notes=[
                "the surface radius is a formula, so this row's surface error is the harness's own; "
                "its useful numbers are the volume/area of the triangulation a solver would receive",
            ],
        )
    )

    # -- mesh local edit: the control --------------------------------------- #
    mesh_started = time.perf_counter()
    base_mesh = _lat_long_surface(lambda direction: DENT_RADIUS_M, rings=rings, segments=segments)
    edited_vertices = []
    moved = 0
    for vertex in base_mesh.vertices:
        distance = math.sqrt(sum(component * component for component in vertex))
        direction = (vertex[0] / distance, vertex[1] / distance, vertex[2] / distance)
        angle = _angle_from_axis(direction)
        if angle <= PROTECTED_FACTOR * DENT_SIGMA_RAD:
            target = dent_radius(direction, DENT_RADIUS_M, DENT_DEPTH_M, DENT_SIGMA_RAD)
            scale = target / distance
            edited_vertices.append((vertex[0] * scale, vertex[1] * scale, vertex[2] * scale))
            moved += 1
        else:
            edited_vertices.append(vertex)
    edited_mesh = Mesh(
        vertices=tuple(edited_vertices),
        faces=base_mesh.faces,
        units=base_mesh.units,
        coordinate_frame=base_mesh.coordinate_frame,
    )
    caster = _RayCaster(edited_mesh)
    mesh_radii: list[float | None] = [caster.radius(direction) for direction in directions]
    base_caster = _RayCaster(base_mesh)
    base_radii: list[float | None] = [base_caster.radius(direction) for direction in directions]
    mesh_seconds = time.perf_counter() - mesh_started
    missing = sum(1 for radius in mesh_radii if radius is None)
    surface_errors = [
        abs(mesh_radii[index] - reference_radii[index])  # type: ignore[operator]
        for index in range(rays)
        if mesh_radii[index] is not None
    ]
    dent_errors = [
        abs(mesh_radii[index] - reference_radii[index])  # type: ignore[operator]
        for index in range(rays)
        if dent_mask[index] and mesh_radii[index] is not None
    ]
    # Two different questions, deliberately kept apart. Compared with the *unedited mesh on
    # the same ray* the protected region measures the edit's leak. Compared with the analytic
    # surface it measures the total error a consumer sees -- which for a triangulation is
    # dominated by chord sag, not by the edit. Only the first one is a gate.
    protected_leaks = [
        abs(mesh_radii[index] - base_radii[index])  # type: ignore[operator]
        for index in range(rays)
        if protected_mask[index] and mesh_radii[index] is not None and base_radii[index] is not None
    ]
    protected_total = [
        abs(mesh_radii[index] - reference_radii[index])  # type: ignore[operator]
        for index in range(rays)
        if protected_mask[index] and mesh_radii[index] is not None
    ]
    transition_leaks = [
        abs(mesh_radii[index] - base_radii[index])  # type: ignore[operator]
        for index in range(rays)
        if transition_mask[index] and mesh_radii[index] is not None and base_radii[index] is not None
    ]
    achieved = [
        DENT_RADIUS_M - mesh_radii[index]  # type: ignore[operator]
        for index in range(rays)
        if dent_mask[index] and mesh_radii[index] is not None
    ]
    # The reference dent is a Gaussian, so its depth averaged over the edit cone is not
    # the depth on the axis. Compare like with like: the measured mean against the
    # reference mean over the same rays, and the axis separately.
    reference_depths = [
        DENT_RADIUS_M - reference_radii[index] for index in range(rays) if dent_mask[index]
    ]
    axis_index = min(range(rays), key=lambda index: angles[index])
    axis_measured_depth = DENT_RADIUS_M - float(mesh_radii[axis_index])  # type: ignore[arg-type]
    # The nearest ray to the axis is not *on* the axis, so the reference depth there is
    # smaller than the peak. Both numbers are reported: the peak is the physical depth,
    # the comparison is the interpolation error at that ray.
    axis_reference_depth = DENT_RADIUS_M - reference_radii[axis_index]
    peak_index = max(
        (index for index in range(rays) if dent_mask[index] and mesh_radii[index] is not None),
        key=lambda index: DENT_RADIUS_M - float(mesh_radii[index]),  # type: ignore[arg-type]
    )
    peak_depth = DENT_RADIUS_M - float(mesh_radii[peak_index])  # type: ignore[arg-type]
    mesh_volume = meshlib.total_signed_volume(edited_mesh)
    mesh_metrics = {
        "dent_depth_measured_mean_m": Metric(
            statistics.fmean(achieved),
            "m",
            reference=f"reference mean over the edit cone ({statistics.fmean(reference_depths):.6f} m)",
        ),
        "dent_depth_mean_error_m": Metric(
            abs(statistics.fmean(achieved) - statistics.fmean(reference_depths)),
            "m",
            reference="reference mean over the same rays",
        ).at_most(0.01 * DENT_DEPTH_M),
        "dent_depth_axis_error_m": Metric(
            abs(axis_measured_depth - axis_reference_depth),
            "m",
            reference=(
                f"reference depth {axis_reference_depth:.6f} m on the ray nearest the axis "
                f"(theta = {angles[axis_index]:.4f} rad)"
            ),
        ).at_most(0.01 * DENT_DEPTH_M),
        "dent_depth_peak_m": Metric(
            peak_depth,
            "m",
            reference=(
                f"target {DENT_DEPTH_M} m at theta = 0; measured on the ray at "
                f"theta = {angles[peak_index]:.4f} rad"
            ),
        ),
        "dent_surface_error_mean_m": Metric(
            statistics.fmean(dent_errors), "m", reference="analytic dent profile"
        ).at_most(0.02 * DENT_DEPTH_M),
        "dent_surface_error_max_m": Metric(max(dent_errors), "m", reference="analytic dent profile"),
        "protected_displacement_mean_m": Metric(
            statistics.fmean(protected_leaks),
            "m",
            reference=(
                "same ray through the unedited mesh: the edit's leak, not the triangulation's sag; "
                f"rays at theta >= {protected_start:.4f} rad"
            ),
        ).at_most(0.02 * DENT_DEPTH_M),
        "protected_displacement_max_m": Metric(
            max(protected_leaks), "m", reference="same ray through the unedited mesh"
        ),
        "protected_total_error_mean_m": Metric(
            statistics.fmean(protected_total),
            "m",
            reference="analytic surface (total error: leak + chord sag, reported not gated)",
        ),
        "surface_error_mean_m": Metric(
            statistics.fmean(surface_errors), "m", reference="analytic surface over every ray"
        ),
        "surface_error_max_m": Metric(
            max(surface_errors), "m", reference="analytic surface over every ray"
        ),
        "transition_band_displacement_m": Metric(
            max(transition_leaks, default=0.0),
            "m",
            reference=(
                f"band {protected_edge:.3f}-{protected_start:.3f} rad where a triangle straddles "
                "the edit support edge (reported, not gated)"
            ),
        ),
        "locality_ratio": Metric(
            max(protected_leaks) / DENT_DEPTH_M,
            "relative",
            reference="worst protected-ray leak / peak edit depth",
        ).at_most(0.05),
        "volume_error": Metric(
            abs(mesh_volume / exact_volume - 1.0), "relative", reference="exact volume"
        ).within(2e-2),
        "discretisation_floor_m": Metric(
            max(
                (
                    DENT_RADIUS_M - float(base_radii[index])  # type: ignore[arg-type]
                    for index in range(rays)
                    if protected_mask[index] and base_radii[index] is not None
                ),
                default=0.0,
            ),
            "m",
            reference=(
                "same grid, unedited sphere: the chord sag a ray between vertices lands on. "
                "The depth gate is 6.0e-4 m, so this is the headroom the gate has"
            ),
        ),
        "ray_coverage": Metric(
            (rays - missing) / rays, "fraction", reference="every ray hits the surface"
        ).at_least(0.999),
        "watertight": Metric(
            1.0 if meshlib.is_watertight(edited_mesh) else 0.0, "bool", reference="closed surface"
        ).at_least(1.0),
        "euler_characteristic": Metric(
            float(
                len(edited_mesh.vertices)
                - len(meshlib.edge_face_counts(edited_mesh))
                + len(edited_mesh.faces)
            ),
            "count",
            reference="2 for a genus-0 closed surface",
        ).equals(2.0),
    }
    representations.append(
        RepresentationResult(
            identifier="mesh_local_edit",
            role="control: explicit triangles, vertices moved in place",
            status=MEASURED,
            metrics=mesh_metrics,
            seconds=mesh_seconds,
            artifacts={"edited_surface": _digest(edited_mesh)},
            notes=[
                f"{moved} of {len(edited_mesh.vertices)} vertices moved; every other vertex is "
                "bit-identical, so the protected-region leak is exactly zero wherever no moved "
                "vertex is reachable, and the edit region pays for linear interpolation instead",
                "no re-projection, no fairing, no smoothing: this is what a plain mesh edit costs",
                "protected_displacement_* is measured against the unedited mesh on the same ray; "
                "protected_total_error_* is measured against the analytic surface and includes "
                "the chord sag of the triangulation",
            ],
        )
    )

    # -- neural field: measured against the same rays, or imported ---------- #
    representations.append(
        _neural_row(
            neural_bundle,
            directions=directions,
            angles=angles,
            dent_mask=dent_mask,
            protected_mask=protected_mask,
            transition_mask=transition_mask,
            reference_radii=reference_radii,
            protected_start=protected_start,
        )
    )

    # -- OCP/STEP: optional kernel path ------------------------------------- #
    representations.append(_ocp_row(ocp_artifact))

    gates = (
        Gate(
            "analytic_row_is_the_reference",
            bool(analytic_metrics["reference_radius_error_max_m"].passed),
            "the harness reproduces the closed form to 1e-12 m",
        ),
        Gate(
            "mesh_edit_hits_the_target_depth",
            bool(mesh_metrics["dent_depth_axis_error_m"].passed),
            f"axis depth error {mesh_metrics['dent_depth_axis_error_m'].value:.3e} m "
            f"(mean over the cone {mesh_metrics['dent_depth_measured_mean_m'].value:.6f} m)",
        ),
        Gate(
            "mesh_edit_is_local",
            bool(mesh_metrics["locality_ratio"].passed),
            f"locality ratio {mesh_metrics['locality_ratio'].value:.3e}",
        ),
        Gate(
            "every_ray_hits_every_surface",
            bool(mesh_metrics["ray_coverage"].passed),
            f"{mesh_metrics['ray_coverage'].value:.4f} coverage",
        ),
        Gate(
            "the_depth_threshold_is_the_same_for_every_measured_row",
            all(
                bool(row.metrics["dent_depth_mean_error_m"].passed)
                for row in representations
                if row.status == MEASURED and "dent_depth_mean_error_m" in row.metrics
            ),
            "1 % of the edit depth, applied unchanged to every representation that measured it "
            "(the analytic row is the reference and is not a candidate)",
        ),
        Gate(
            "no_representation_went_unmeasured_silently",
            all(
                result.status in (MEASURED, IMPORTED, NOT_SUPPORTED)
                for result in representations
            ),
            "every row is measured, imported or explicitly not_supported",
        ),
    )
    return CaseResult(
        identifier="case_a_sphere_dent",
        question=(
            "A local edit with a closed-form reference: how much of the error is the edit, "
            "and how much is the representation?"
        ),
        reference=(
            f"sphere R = {DENT_RADIUS_M} m with a Gaussian dent of depth {DENT_DEPTH_M} m and "
            f"sigma {DENT_SIGMA_RAD} rad on +x; {rays} stratified rays, dent < 1 sigma, "
            f"protected > {PROTECTED_FACTOR} sigma + 1.5 cells"
        ),
        representations=tuple(representations),
        gates=gates,
        notes=[
            "the mesh row moves only vertices inside 2.5 sigma and does not re-project: the error "
            "it shows is interpolation between vertices, not a defect of meshes",
            "every row is probed along the same stratified ray set: a third of the rays lie inside "
            "the edit cone so a small edit is resolved, the rest are uniform, so no row sees an "
            "easier sample",
            "an unweighted mean over rays is not a volume average: the rows are compared with each "
            "other, not with a surface integral",
        ],
        seconds=time.perf_counter() - started,
    )


@lru_cache(maxsize=1)
def have_neural() -> bool:
    """True when the optional neural runtime can actually be imported here.

    `find_spec` is not enough: a partially installed stack (torch present, safetensors
    missing) finds a spec and then fails at import, which is exactly the state this project
    was in once. The honest question is whether the import succeeds.
    """
    for name in ("numpy", "torch", "safetensors"):
        try:
            importlib.import_module(name)
        except ImportError:
            return False
    return True


def _field_radii(
    evaluate: Callable[[Any], Any],
    directions: Sequence[Sequence[float]],
    *,
    r_min: float = 0.30,
    r_max: float = 0.70,
    samples: int = 241,
    iterations: int = 48,
) -> list[float | None]:
    """Where a truncated signed distance field crosses zero along each ray.

    A learned field has no triangles to intersect, so it has to be marched: sample the ray,
    take the first inside-to-outside crossing, then bisect it to float precision. All rays
    are bisected together (one batched call per iteration) because the network is the
    expensive part and one point at a time would take minutes.
    """
    import numpy as np

    rays = np.asarray(directions, dtype=np.float64).reshape(-1, 3)
    grid = np.linspace(r_min, r_max, samples)
    points = (rays[:, None, :] * grid[None, :, None]).reshape(-1, 3)
    values = np.asarray(evaluate(points), dtype=np.float64).reshape(len(rays), samples)
    low = np.full(len(rays), np.nan)
    high = np.full(len(rays), np.nan)
    for index, row in enumerate(values):
        crossing = np.nonzero((row[:-1] <= 0.0) & (row[1:] >= 0.0))[0]
        if crossing.size:
            low[index] = grid[int(crossing[0])]
            high[index] = grid[int(crossing[0]) + 1]
    found = np.nonzero(~np.isnan(low))[0]
    for _ in range(iterations):
        middle = 0.5 * (low[found] + high[found])
        mid_values = np.asarray(evaluate(rays[found] * middle[:, None]), dtype=np.float64)
        inside = mid_values <= 0.0
        low[found[inside]] = middle[inside]
        high[found[~inside]] = middle[~inside]
    hit = 0.5 * (low + high)
    return [None if np.isnan(value) else float(value) for value in hit]


def _neural_row(
    bundle: str | Path | None,
    *,
    directions: Sequence[Sequence[float]],
    angles: Sequence[float],
    dent_mask: Sequence[bool],
    protected_mask: Sequence[bool],
    transition_mask: Sequence[bool],
    reference_radii: Sequence[float],
    protected_start: float,
) -> RepresentationResult:
    """A verified ``.ngo`` bundle, queried on exactly the rays the other rows were.

    Two kinds of number appear here and they are labelled. *Measured* numbers come from
    this process asking the field where the surface is, on the same rays and against the
    same closed form as every other row. *Imported* numbers come from the run that produced
    the bundle, because the bundle's own validation compared against its parent field, which
    is not shipped -- so the leak of the edit is imported, and clearly marked as such.
    """
    identifier = "neural_field"
    role = "learned representation: versioned .ngo object, queried as an implicit field"
    if bundle is None:
        zero = "no bundle supplied; produce one with `python examples/neural_geometry/edit_sphere_dent.py`"
        return RepresentationResult(
            identifier=identifier,
            role=role,
            status=NOT_SUPPORTED,
            metrics={
                "dent_surface_error_mean_m": _missing(NOT_SUPPORTED, zero, "m"),
                "protected_displacement_mean_m": _missing(NOT_SUPPORTED, zero, "m"),
                "protected_total_error_mean_m": _missing(NOT_SUPPORTED, zero, "m"),
                "locality_ratio": _missing(NOT_SUPPORTED, zero, "relative"),
            },
            notes=[
                "this row is deliberately not simulated: a learned field's numbers must come from a "
                "trained artifact, and an untrained stand-in would be a mock result",
            ],
        )
    path = Path(bundle)
    from .neural_geometry import load_bundle

    try:
        manifest = load_bundle(path)
    except Exception as exc:  # a benchmark records a failed artifact instead of raising
        return RepresentationResult(
            identifier=identifier,
            role=role,
            status=FAILED,
            metrics={
                "dent_surface_error_mean_m": Metric(
                    value=None, unit="m", status=FAILED, note=f"bundle did not load: {exc}"
                )
            },
            notes=[f"{path} failed verification; nothing from it was used"],
        )
    validation = json.loads((path / manifest.validation_file).read_text(encoding="utf-8"))
    reported = validation.get("metrics", {})
    artifacts = {"bundle": str(path), "weights_sha256": manifest.weights.sha256}
    imported: dict[str, Metric] = {}
    # Two producers write bundles today (`fit_sphere`/`apply_dent_edit` and the
    # edit_sphere_dent example) and they name the dent error differently. Read either,
    # and say which key was used rather than guessing silently.
    reported_dent_key = next(
        (key for key in ("dent_surface_error_mean_m", "dent_region_error_mean_m") if key in reported),
        None,
    )
    if reported_dent_key is not None:
        imported["reported_dent_surface_error_mean_m"] = Metric(
            float(reported[reported_dent_key]),
            "m",
            status=IMPORTED,
            reference="the producing run's own validation, on its own ray set",
            note=(
                f"{reported_dent_key} from {manifest.validation_file}; compare with this "
                "harness's measurement"
            ),
        )
    if "protected_mean_displacement_m" in reported:
        imported["protected_displacement_mean_m"] = Metric(
            float(reported["protected_mean_displacement_m"]),
            "m",
            status=IMPORTED,
            reference="displacement against the producing run's parent field (the edit's leak)",
            note="the parent weights are not part of the bundle, so this cannot be re-measured here",
        )
        imported["locality_ratio"] = Metric(
            float(reported["protected_mean_displacement_m"]) / DENT_DEPTH_M,
            "relative",
            status=IMPORTED,
            reference="imported leak / peak edit depth",
        )
    producing_failed = validation.get("status") != "passed"

    if not have_neural():
        return RepresentationResult(
            identifier=identifier,
            role=role,
            status=IMPORTED,
            metrics=imported,
            artifacts=artifacts,
            notes=[
                "imported, not re-measured: the neural extra is not installed in this process "
                '(`pip install "tensorcad[neural]"`), so the numbers come from the run that produced '
                "the bundle; its hashes were verified here",
                "with the extra installed this row is re-measured on this case's own rays instead",
                *(["the producing run reports status != passed"] if producing_failed else []),
            ],
        )

    started = time.perf_counter()
    try:
        from .neural.objects import field_from_object

        field = field_from_object(path)
        radii = _field_radii(field.numpy_callable(), directions)
    except Exception as exc:  # a broken artifact is a failed row, not a crashed benchmark
        return RepresentationResult(
            identifier=identifier,
            role=role,
            status=FAILED,
            metrics={
                "dent_surface_error_mean_m": Metric(
                    value=None, unit="m", status=FAILED, note=f"field could not be queried: {exc}"
                )
            },
            artifacts=artifacts,
            notes=[f"{path} verified but is not evaluable as a field: {exc}"],
        )
    seconds = time.perf_counter() - started

    rays = len(directions)
    missing = sum(1 for radius in radii if radius is None)
    errors = [
        abs(radii[index] - reference_radii[index])  # type: ignore[operator]
        for index in range(rays)
        if radii[index] is not None
    ]
    dent_errors = [
        abs(radii[index] - reference_radii[index])  # type: ignore[operator]
        for index in range(rays)
        if dent_mask[index] and radii[index] is not None
    ]
    protected_errors = [
        abs(radii[index] - reference_radii[index])  # type: ignore[operator]
        for index in range(rays)
        if protected_mask[index] and radii[index] is not None
    ]
    transition_errors = [
        abs(radii[index] - reference_radii[index])  # type: ignore[operator]
        for index in range(rays)
        if transition_mask[index] and radii[index] is not None
    ]
    achieved = [
        DENT_RADIUS_M - radii[index]  # type: ignore[operator]
        for index in range(rays)
        if dent_mask[index] and radii[index] is not None
    ]
    reference_depths = [
        DENT_RADIUS_M - reference_radii[index] for index in range(rays) if dent_mask[index]
    ]
    axis_index = min(range(rays), key=lambda index: angles[index])
    axis_measured_depth = DENT_RADIUS_M - float(radii[axis_index])  # type: ignore[arg-type]
    axis_reference_depth = DENT_RADIUS_M - reference_radii[axis_index]
    peak_index = max(
        (index for index in range(rays) if dent_mask[index] and radii[index] is not None),
        key=lambda index: DENT_RADIUS_M - float(radii[index]),  # type: ignore[arg-type]
    )
    measured = {
        "dent_depth_measured_mean_m": Metric(
            statistics.fmean(achieved),
            "m",
            reference=f"reference mean over the edit cone ({statistics.fmean(reference_depths):.6f} m)",
        ),
        "dent_depth_mean_error_m": Metric(
            abs(statistics.fmean(achieved) - statistics.fmean(reference_depths)),
            "m",
            reference="reference mean over the same rays",
        ).at_most(0.01 * DENT_DEPTH_M),
        "dent_depth_axis_error_m": Metric(
            abs(axis_measured_depth - axis_reference_depth),
            "m",
            reference=(
                f"reference depth {axis_reference_depth:.6f} m on the ray nearest the axis "
                f"(theta = {angles[axis_index]:.4f} rad)"
            ),
        ),
        "dent_depth_peak_m": Metric(
            DENT_RADIUS_M - float(radii[peak_index]),  # type: ignore[arg-type]
            "m",
            reference=f"target {DENT_DEPTH_M} m on the axis; measured on the ray at theta = "
            f"{angles[peak_index]:.4f} rad",
        ),
        "dent_surface_error_mean_m": Metric(
            statistics.fmean(dent_errors), "m", reference="analytic dent profile"
        ),
        "dent_surface_error_max_m": Metric(max(dent_errors), "m", reference="analytic dent profile"),
        "protected_total_error_mean_m": Metric(
            statistics.fmean(protected_errors),
            "m",
            reference=(
                "analytic surface outside the edit cone: for a learned field this is the fit error "
                "it carries everywhere, not the edit's leak (reported, not gated)"
            ),
        ),
        "protected_total_error_max_m": Metric(max(protected_errors), "m", reference="analytic surface"),
        "transition_band_total_error_m": Metric(
            max(transition_errors, default=0.0),
            "m",
            reference="band between the edit and the protected region",
        ),
        "surface_error_mean_m": Metric(
            statistics.fmean(errors), "m", reference="analytic surface over every ray"
        ),
        "surface_error_max_m": Metric(max(errors), "m", reference="analytic surface over every ray"),
        "ray_coverage": Metric(
            (rays - missing) / rays, "fraction", reference="every ray crosses the surface"
        ).at_least(0.999),
        "truncation_distance_m": Metric(
            float(field.truncation_distance), "m", reference="the field saturates outside this band"
        ),
    }
    metrics = {**measured, **imported}
    notes = [
        f"queried on the same {rays} rays as every other row, by marching the field and "
        f"bisecting the zero crossing ({field.parameters()['architecture_id']}, "
        f"truncation {field.truncation_distance} m)",
        "dent_depth_* and dent_surface_error_* are measured here, not imported: the learned "
        "surface is compared with the same closed form as the mesh row",
        "the protected region is only measured as *total* error (fit error plus leak); the "
        "leak alone is the imported number, because the parent weights are not in the bundle",
        "no gate is attached to this row's accuracy: no published tolerance exists for a learned "
        "field, and inventing one after seeing these numbers would be tuning, not measurement",
        *(["the producing run reports status != passed"] if producing_failed else []),
    ]
    return RepresentationResult(
        identifier=identifier,
        role=role,
        status=MEASURED,
        metrics=metrics,
        seconds=seconds,
        artifacts=artifacts,
        notes=notes,
    )


def _ocp_row(artifact: str | Path | None) -> RepresentationResult:
    """A faceted STEP round trip, when OpenCASCADE is installed."""
    from .adapters.ocp import OcpUnavailable, have_ocp, read_step, write_step

    if not have_ocp():
        return RepresentationResult(
            identifier="ocp_step",
            role="kernel path: faceted solid exported as STEP",
            status=NOT_SUPPORTED,
            metrics={
                "step_volume_error": _missing(
                    NOT_SUPPORTED, 'needs the optional wheels: pip install "tensorcad[ocp]"', "relative"
                )
            },
            notes=["a faceted solid is not a feature B-rep; when this row runs, it reports the facet "
                   "count next to the round-trip volume so the distinction is measurable"],
        )
    path = Path(artifact) if artifact is not None else Path("results/bench/sphere_dent.step")
    mesh = _lat_long_surface(
        lambda direction: dent_radius(direction, DENT_RADIUS_M, DENT_DEPTH_M, DENT_SIGMA_RAD),
        rings=20,
        segments=40,
    )
    started = time.perf_counter()
    try:
        record = write_step(mesh, path, round_trip=True)
    except OcpUnavailable as exc:  # pragma: no cover - have_ocp() already checked
        return RepresentationResult(
            identifier="ocp_step", role="kernel path", status=NOT_SUPPORTED,
            metrics={"step_volume_error": _missing(NOT_SUPPORTED, str(exc), "relative")},
        )
    seconds = time.perf_counter() - started
    imported = read_step(path)
    exact_volume = _exact_volume(DENT_RADIUS_M, DENT_DEPTH_M, DENT_SIGMA_RAD)
    volume_error = abs(record.volume_m3 / exact_volume - 1.0)
    return RepresentationResult(
        identifier="ocp_step",
        role="kernel path: faceted solid exported as STEP",
        status=MEASURED,
        metrics={
            "step_volume_error": Metric(volume_error, "relative", reference="exact volume").within(2e-2),
            "round_trip_volume_error": Metric(
                abs(record.round_trip_volume_error or 0.0), "relative", reference="written volume"
            ).within(1e-9),
            "faceted_solids": Metric(float(imported.solids), "count", reference="one closed shell"),
            "planar_faces": Metric(
                float(imported.faces),
                "count",
                reference="one per triangle: evidence that this is a tessellation",
                note=f"{len(mesh.faces)} triangles in, {imported.faces} planar faces out",
            ),
            "file_bytes": Metric(float(record.file_bytes), "bytes", reference="artifact size"),
        },
        seconds=seconds,
        artifacts={"step": str(path)},
        notes=["sewing is not feature modelling: every face in this file is planar"],
    )


# --------------------------------------------------------------------------- #
# Case B: the parametric hull, from target to exported mesh
# --------------------------------------------------------------------------- #


def hull_handoff_case(
    *,
    length: float = 12.0,
    beam: float = 2.4,
    draft: float = 0.9,
    write_artifacts: bool = False,
    artifact_dir: str | Path | None = None,
) -> CaseResult:
    """Case B: a design target, a converged hull, and the artifact a solver receives.

    The reference hull's *analytic* volume is the ground truth (``hull.analytic_volume``),
    so the mesh error and the design-loop residual are measured against something that is
    not the mesh itself. The export round trip is measured by re-reading the gmsh file.
    """
    started = time.perf_counter()
    from .adapters.fem import fem_surface, read_gmsh, write_gmsh
    from .design import DesignProblem, DesignTargets
    from .hull import HullParameters, analytic_volume, build_mesh

    reference = HullParameters(length=length, beam=beam, draft=draft)
    exact = analytic_volume(length, beam, draft)
    reference_mesh, _ = build_mesh(reference)
    measured = meshlib.total_signed_volume(reference_mesh)
    area = meshlib.surface_area(reference_mesh)

    design_started = time.perf_counter()
    solution = DesignProblem(
        base=HullParameters(length=0.85 * length, beam=0.8 * beam, draft=0.9 * draft),
        targets=DesignTargets(
            displacement_m3=exact,
            beam_to_draft=beam / draft,
            length_to_beam=length / beam,
        ),
        variables=("length", "beam", "draft"),
    ).solve()
    design_seconds = time.perf_counter() - design_started

    export_started = time.perf_counter()
    surface = fem_surface(reference_mesh)
    # The round trip needs a real file on disk. When the caller did not ask for artifacts
    # it goes to a scratch directory, so publishing a result does not scatter .msh files
    # next to the report.
    persistent = Path(artifact_dir) if artifact_dir is not None else None
    with tempfile.TemporaryDirectory() as scratch:
        target = persistent if persistent is not None else Path(scratch)
        gmsh_path = write_gmsh(surface, target / "hull_surface.msh")
        parsed = read_gmsh(gmsh_path)
    export_seconds = time.perf_counter() - export_started
    worst_node = max(
        abs(a - b)
        for node_id, node in enumerate(surface.nodes, start=1)
        for a, b in zip(parsed["node_table"][node_id], node, strict=True)
    )
    area_round_trip = abs(parsed["area_by_region_m2"][1] - area) / area

    metrics = {
        "volume_error": Metric(
            abs(measured / exact - 1.0), "relative", reference="analytic_volume"
        ).within(1e-2),
        "design_loop_volume_error": Metric(
            abs(solution.final.hydrostatics.displacement_volume / exact - 1.0),
            "relative",
            reference="analytic_volume",
        ).within(1e-6),
        "design_loop_solved": Metric(
            1.0 if solution.solved else 0.0, "bool", reference="converged and verified"
        ).at_least(1.0),
        "design_loop_evaluations": Metric(
            float(len(solution.records)), "count", reference="iteration budget"
        ),
        "design_loop_seconds": Metric(design_seconds, "s", reference="wall clock around solve()"),
        "watertight": Metric(1.0 if meshlib.is_watertight(reference_mesh) else 0.0, "bool").at_least(1.0),
        "euler_characteristic": Metric(
            float(
                len(reference_mesh.vertices)
                - len(meshlib.edge_face_counts(reference_mesh))
                + len(reference_mesh.faces)
            ),
            "count",
            reference="2",
        ).at_most(2.0),
        "export_node_error_m": Metric(worst_node, "m", reference="exact round trip").at_most(0.0),
        "export_area_error": Metric(area_round_trip, "relative", reference="mesh area").within(1e-12),
        "surface_area_m2": Metric(area, "m^2", reference="integrated from the mesh"),
    }

    # Grid dependence, measured rather than assumed: same parameters, twice the density.
    refined = HullParameters(
        length=length, beam=beam, draft=draft, nx=2 * reference.nx - 1, nz=2 * reference.nz - 1
    )
    refined_mesh, _ = build_mesh(refined)
    refined_volume = meshlib.total_signed_volume(refined_mesh)
    metrics["grid_refinement_volume_change"] = Metric(
        abs(refined_volume / measured - 1.0), "relative", reference=f"({reference.nx}, {reference.nz})"
    ).within(5e-2)

    artifacts = {"hull_surface": _digest(reference_mesh)}
    if write_artifacts and persistent is not None:
        from .export import write_obj

        obj_path = write_obj(reference_mesh, persistent / "hull_surface.obj")
        artifacts["hull_surface_obj_sha256"] = hashlib.sha256(obj_path.read_bytes()).hexdigest()
        artifacts["gmsh"] = str(gmsh_path)
    representations = (
        RepresentationResult(
            identifier="parametric_hull",
            role="engineering case: three design variables, analytic volume reference",
            status=MEASURED,
            metrics=metrics,
            seconds=design_seconds + export_seconds,
            artifacts=artifacts,
            notes=[
                f"the design loop converged from 85 % of the target dimensions in "
                f"{len(solution.records)} evaluations ({solution.status})",
                f"analytic volume {exact:.6f} m^3 vs mesh {measured:.6f} m^3 at "
                f"({reference.nx}, {reference.nz}); refining to ({refined.nx}, {refined.nz}) moves it "
                f"by {100 * abs(refined_volume / measured - 1.0):.3f} %",
            ],
        ),
        _ocp_hull_row(reference_mesh, exact, persistent if write_artifacts else None),
    )
    gates = (
        Gate(
            "hull_mesh_matches_the_analytic_volume",
            bool(metrics["volume_error"].passed),
            f"{100 * (metrics['volume_error'].value or 0.0):.3f} % error at "
            f"({reference.nx}, {reference.nz})",
        ),
        Gate(
            "design_loop_reaches_the_target",
            bool(metrics["design_loop_volume_error"].passed) and bool(metrics["design_loop_solved"].passed),
            f"{solution.status}, {len(solution.records)} evaluations",
        ),
        Gate(
            "export_round_trip_is_exact",
            bool(metrics["export_node_error_m"].passed) and bool(metrics["export_area_error"].passed),
            f"worst node displacement {worst_node:.1e} m, area {area_round_trip:.1e}",
        ),
    )
    return CaseResult(
        identifier="case_b_hull_handoff",
        question=(
            "The engineering case: can a target become a verified hull, and does the artifact a "
            "solver reads describe the same surface?"
        ),
        reference=f"L = {length} m, B = {beam} m, T = {draft} m; analytic volume {exact:.9f} m^3",
        representations=representations,
        gates=gates,
        notes=[
            "the mesh volume here is a discretisation of an analytic form, so its error is the grid's: "
            "the refinement row measures how much",
        ],
        seconds=time.perf_counter() - started,
    )


def _ocp_hull_row(mesh: Mesh, exact_volume: float, target: Path | None) -> RepresentationResult:
    from .adapters.ocp import have_ocp

    if not have_ocp():
        return RepresentationResult(
            identifier="ocp_step",
            role="kernel path: faceted solid exported as STEP",
            status=NOT_SUPPORTED,
            metrics={
                "step_volume_error": _missing(
                    NOT_SUPPORTED, 'needs the optional wheels: pip install "tensorcad[ocp]"', "relative"
                )
            },
        )
    from .adapters.ocp import read_step, write_step

    path = (target or Path("results/bench")) / "hull_solid.step"
    started = time.perf_counter()
    record = write_step(mesh, path, round_trip=True)
    imported = read_step(path)
    seconds = time.perf_counter() - started
    return RepresentationResult(
        identifier="ocp_step",
        role="kernel path: faceted solid exported as STEP",
        status=MEASURED,
        metrics={
            "step_volume_error": Metric(
                abs(record.volume_m3 / exact_volume - 1.0), "relative", reference="analytic volume"
            ).within(1e-2),
            "round_trip_volume_error": Metric(
                abs(record.round_trip_volume_error or 0.0), "relative", reference="written volume"
            ).within(1e-9),
            "planar_faces": Metric(float(imported.faces), "count", reference="one per triangle"),
            "file_bytes": Metric(float(record.file_bytes), "bytes"),
        },
        seconds=seconds,
        artifacts={"step": str(path)},
    )


# --------------------------------------------------------------------------- #
# the runner
# --------------------------------------------------------------------------- #


def run_benchmark(
    *,
    quick: bool = False,
    neural_bundle: str | Path | None = None,
    ocp_artifact: str | Path | None = None,
    write_artifacts: bool = False,
    artifact_dir: str | Path | None = None,
    command: str = "",
) -> BenchmarkResult:
    """Run both cases under one environment record and one verdict.

    ``quick`` shrinks the sampling (fewer rays, coarser rings) for CI. It does not relax a
    gate: a quick run that fails the depth or locality check is a failed run, because the
    gate is what protects the number, not the sample size.
    """
    started = time.perf_counter()
    rays = 128 if quick else 512
    # `quick` reduces the *sampling*, never the surface: the grid sets the chord-sag floor
    # that the depth gate is compared against, so a coarser quick grid would move the thing
    # being measured (the ladder is in docs/benchmarks/representation-comparison.md). At
    # 80x160 the floor is ~1.4e-4 m against a 6.0e-4 m gate, a factor of four in hand.
    rings, segments = 80, 160
    cases = (
        sphere_dent_case(
            rays=rays,
            rings=rings,
            segments=segments,
            neural_bundle=neural_bundle,
            ocp_artifact=ocp_artifact,
        ),
        hull_handoff_case(write_artifacts=write_artifacts, artifact_dir=artifact_dir),
    )
    return BenchmarkResult(
        cases=cases,
        environment=_environment(),
        seconds=time.perf_counter() - started,
        command=command,
    )
