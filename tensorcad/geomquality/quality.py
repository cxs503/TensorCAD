"""Canonical element-quality metrics for triangle meshes, and the gate that uses them.

Why this module exists.  Before it, this repository answered "how good is this
triangulation?" in four mutually incompatible dialects, all of them called *aspect
ratio*:

==========================================  ============================  =============
dialect                                     quantity actually computed    equilateral
==========================================  ============================  =============
:func:`tensorcad.gate._facet_aspect`        ``l_max**2 / A``               2.3094
:func:`tensorcad.validation.measure_quality` ``l_max**2 / (2A)``           1.1547
:func:`tensorcad.parametric_hull.mesh_quality_report`  ``l_max / l_min``    1.0000
:mod:`tensorcad.mesh_conditioning`          ``l_max / l_min`` (guard 25)   1.0000
==========================================  ============================  =============

Three of the four are *not* an aspect ratio in the sense any mesher, any solver
preprocessor or any of the upstream libraries a shipyard already owns will report, and
two of them differ from each other by exactly a factor of two.  A number that cannot be
compared with a number produced by another tool is not a quality metric; it is an
in-house convention with a misleading name.  ``docs/mesh-quality-lib.md`` measures the
four dialects against this module on the same three meshes, so the discrepancy is a
number a reader can check rather than a claim.

The metric set below is the one :func:`meshmagick.Mesh.print_quality` reports --
meshmagick delegates to ``vtk.vtkMeshQuality``, which is a transcription of the Verdict
library -- plus the two Gmsh measures a CAD/CFD handoff actually argue about.  Every
definition is written out in this file as a formula, every one is normalised so that the
*ideal element* scores 1 (except the condition number, whose ideal is ``1/sqrt(3)`` and
which is documented as such), and every one is invariant under rotation, translation
and uniform scaling, so the value is comparable across hull sizes.

Upstream-provenance note, deliberately blunt.  This environment has no VTK, no Gmsh
binary and no meshmagick *executable* (meshmagick's own quality path needs VTK; the
`sources <https://gitlab.com/Novanum/MeshMagick>`_ were read, and
``meshmagick/mesh.py:459-540`` is the wiring).  The formulas below are therefore
transcribed, not executed against the upstream tools.  What is *not* a transcription:

* every measure is checked against a closed-form value on canonical elements
  (equilateral, right isosceles, 3-4-5, degenerate) to ``1e-12``;
* every measure is recomputed by an independent second path -- exact rational
  arithmetic on rational meshes, ``numpy.linalg.svd`` instead of the algebraic
  identity for the conditioning, the circumradius/inradius route instead of the
  edge-length route for the aspect ratio -- and the two paths must agree to ``1e-6``
  on the hull meshes (``tools/mesh_quality_crosscheck.py``);
* ``tools/quality_upstream_oracle.py`` runs the real thing and reports the difference
  the moment ``vtk`` or ``gmsh`` is importable, so a host that has them turns the
  transcription into a measurement without touching this file.

:data:`TRIANGLE_METRICS` is the whole vocabulary; :func:`mesh_quality` is the aggregate;
:func:`quality_gate` is the CI door.  Quadrilateral measures are *not* implemented and
this module says so instead of guessing: :class:`tensorcad.geometry.Mesh` only stores
triangles (``Mesh.__post_init__`` raises on any other arity), so a quad metric here
would be dead code with an untested definition.
"""
from __future__ import annotations

import json
import math
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any

from ..geometry import Mesh

Vec3 = tuple[float, float, float]

__all__ = [
    "DEFAULT_CRITERIA",
    "IDEAL_VALUES",
    "METRIC_COMPARISON",
    "METRIC_DEFINITIONS",
    "QualityCheck",
    "QualityVerdict",
    "TriangleQuality",
    "mesh_quality",
    "quality_gate",
    "triangle_quality",
]

_SQRT3 = math.sqrt(3.0)
_FOUR_SQRT3 = 4.0 * _SQRT3

#: Comparison sense of each measure: ``"min"`` means "the smallest element value must be
#: at least the limit", ``"max"`` means "the largest element value must be at most it".
METRIC_COMPARISON: dict[str, str] = {
    "min_angle_deg": "min",
    "max_angle_deg": "max",
    "gamma": "min",
    "aspect_ratio": "max",
    "aspect_frobenius": "max",
    "edge_ratio": "max",
    "radius_ratio": "min",
    "scaled_jacobian": "min",
    "condition_number_inverse": "min",
    "area": "min",
}

#: Value of each measure on the ideal (equilateral) triangle.  ``inf``/``0`` where the
#: ideal value is not a useful reference.
IDEAL_VALUES: dict[str, float] = {
    "min_angle_deg": 60.0,
    "max_angle_deg": 60.0,
    "gamma": 1.0,
    "aspect_ratio": 1.0,
    "aspect_frobenius": 1.0,
    "edge_ratio": 1.0,
    "radius_ratio": 1.0,
    "scaled_jacobian": _SQRT3 / 2.0,
    "condition_number_inverse": 1.0 / _SQRT3,
}

#: One line per measure: what it is, in what unit, and where the definition comes from.
METRIC_DEFINITIONS: dict[str, str] = {
    "area": "A, the triangle area (model units squared).",
    "min_angle_deg": "smallest interior angle, degrees.",
    "max_angle_deg": "largest interior angle, degrees.",
    "edge_ratio": "l_max/l_min -- Verdict 'Edge Ratio' (tri), ideal 1.",
    "aspect_ratio": (
        "l_max*(l_a+l_b+l_c)/(4*sqrt(3)*A) -- the edge-length form of Verdict 'Aspect "
        "Ratio' (tri), ideal 1. NOT the radius form R/(2r): the two are not equal for a "
        "triangle (3-4-5 example in docs/mesh-quality-lib.md: 1.4434 vs 1.2500), so the "
        "radius form is reported separately as 1/radius_ratio."
    ),
    "aspect_frobenius": (
        "(l_a^2+l_b^2+l_c^2)/(4*sqrt(3)*A) -- Verdict 'Aspect Frobenius' (tri), ideal 1."
    ),
    "gamma": "4*sqrt(3)*A/(l_a^2+l_b^2+l_c^2) -- Gmsh 'gamma' (quality), ideal 1.",
    "radius_ratio": (
        "2*r/R = 8*A^2/(s*l_a*l_b*l_c) with s the semiperimeter -- Verdict 'Radius "
        "Ratio' (tri) layout, ideal 1."
    ),
    "scaled_jacobian": (
        "2*A/(l_mid*l_max) -- scaled Jacobian of the reference->physical map "
        "(a.k.a. min scaled Jacobian), ideal sqrt(3)/2."
    ),
    "condition_number_inverse": (
        "sigma_min/sigma_max of the 2x2 Jacobian [[v1-v0],[v2-v0]] (Gmsh 'minSICN' "
        "family), computed from the two invariants sigma1*sigma2 = 2A and "
        "sigma1^2+sigma2^2 = |u|^2+|v|^2; ideal 1/sqrt(3)."
    ),
}

#: Default gate, and why each limit is what it is.  The two angle/quality limits are the
#: thresholds this repository already ships elsewhere (``validation.measure_quality``
#: calls a triangle thinner than 2 degrees a sliver; ``mesh_conditioning`` refuses a mesh
#: above aspect ratio 25).  Both are re-expressed here in the canonical definitions, and
#: ``docs/mesh-quality-lib.md`` reports the margin on clean hulls versus injected slivers
#: so the numbers can be argued with.
DEFAULT_CRITERIA: dict[str, float] = {
    "min_angle_deg": 2.0,
    "gamma": 0.05,
    "aspect_ratio": 25.0,
    "edge_ratio": 100.0,
}

#: What to do about a failing measure, in the vocabulary of :mod:`.defects`.
_REPAIR_HINTS: dict[str, str] = {
    "min_angle_deg": "smooth or collapse the cap/sliver at the reported element",
    "gamma": "remesh locally; a gamma below 0.05 is a cap, not a facet",
    "aspect_ratio": "split the long edge or collapse the short one",
    "edge_ratio": "collapse the short edge (normally a micro edge at a sharp feature)",
    "radius_ratio": "equilateralise the element (Laplacian/Taubin + reconnection)",
    "scaled_jacobian": "untangle the element before compression",
    "condition_number_inverse": "improve element conditioning (local reconnection)",
    "area": "remove zero-area facets before the solver sees them",
}


@dataclass(frozen=True)
class TriangleQuality:
    """Every measure of :data:`METRIC_DEFINITIONS` for one triangle.

    ``degenerate`` is true when the area is zero at working precision; every ratio is
    then reported as ``inf``/``0`` in the direction that makes the element the worst
    possible one, never as ``nan`` -- a ``nan`` silently passes a ``min``/``max`` gate.
    """

    index: int
    area: float
    edge_lengths: tuple[float, float, float]
    min_angle_deg: float
    max_angle_deg: float
    edge_ratio: float
    aspect_ratio: float
    aspect_frobenius: float
    gamma: float
    radius_ratio: float
    scaled_jacobian: float
    condition_number_inverse: float
    degenerate: bool

    def value(self, metric: str) -> float:
        try:
            return float(getattr(self, metric))
        except AttributeError as error:  # pragma: no cover - guarded by METRIC_DEFINITIONS
            raise KeyError(f"unknown metric {metric!r}") from error

    def to_dict(self) -> dict[str, Any]:
        data: dict[str, Any] = {
            "index": self.index,
            "area": self.area,
            "edge_lengths": list(self.edge_lengths),
            "degenerate": self.degenerate,
        }
        data.update({name: self.value(name) for name in METRIC_DEFINITIONS})
        return data


def _angle_at(apex: Vec3, first: Vec3, second: Vec3) -> float:
    """Interior angle at ``apex``, in degrees, from the two unit edge directions."""
    ux, uy, uz = first[0] - apex[0], first[1] - apex[1], first[2] - apex[2]
    vx, vy, vz = second[0] - apex[0], second[1] - apex[1], second[2] - apex[2]
    nu = math.sqrt(ux * ux + uy * uy + uz * uz)
    nv = math.sqrt(vx * vx + vy * vy + vz * vz)
    if nu <= 0.0 or nv <= 0.0:
        return 0.0
    cosine = (ux * vx + uy * vy + uz * vz) / (nu * nv)
    return math.degrees(math.acos(max(-1.0, min(1.0, cosine))))


def triangle_quality(
    a: Vec3,
    b: Vec3,
    c: Vec3,
    *,
    index: int = 0,
) -> TriangleQuality:
    """All measures of the triangle ``(a, b, c)``.

    The formulas are the ones in :data:`METRIC_DEFINITIONS`; nothing is normalised by
    a reference element of the same *size*, so the result is scale free apart from
    ``area``.
    """
    lab = math.dist(a, b)
    lbc = math.dist(b, c)
    lca = math.dist(c, a)
    ux, uy, uz = b[0] - a[0], b[1] - a[1], b[2] - a[2]
    vx, vy, vz = c[0] - a[0], c[1] - a[1], c[2] - a[2]
    nx, ny, nz = uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx
    two_area = math.sqrt(nx * nx + ny * ny + nz * nz)
    area = 0.5 * two_area
    scale = max(lab, lbc, lca)
    # A triangle whose area is below round-off relative to its longest edge has no
    # meaningful shape measures; report it as the worst element instead of dividing by it.
    degenerate = area <= 0.0 or area <= 1e-15 * scale * scale
    if degenerate:
        return TriangleQuality(
            index=index,
            area=area,
            edge_lengths=(lab, lbc, lca),
            min_angle_deg=0.0,
            max_angle_deg=180.0,
            edge_ratio=math.inf,
            aspect_ratio=math.inf,
            aspect_frobenius=math.inf,
            gamma=0.0,
            radius_ratio=0.0,
            scaled_jacobian=0.0,
            condition_number_inverse=0.0,
            degenerate=True,
        )

    short, mid, long_ = sorted((lab, lbc, lca))
    perimeter = lab + lbc + lca
    squares = lab * lab + lbc * lbc + lca * lca
    semiperimeter = 0.5 * perimeter
    frobenius = squares / (_FOUR_SQRT3 * area)
    # |u|^2 + |v|^2 is the sum of the squared singular values of [[u],[v]] and 2A is
    # their product, so the condition number needs no SVD here (verified against one in
    # tools/mesh_quality_crosscheck.py).
    trace = ux * ux + uy * uy + uz * uz + vx * vx + vy * vy + vz * vz
    product = two_area
    discriminant = max(0.0, trace * trace - 4.0 * product * product)
    sigma_max_sq = 0.5 * (trace + math.sqrt(discriminant))
    sigma_min_sq = 0.5 * (trace - math.sqrt(discriminant))
    sicn = math.sqrt(sigma_min_sq / sigma_max_sq) if sigma_max_sq > 0.0 else 0.0
    radius_ratio = (
        8.0 * area * area / (semiperimeter * lab * lbc * lca) if semiperimeter > 0.0 else 0.0
    )
    return TriangleQuality(
        index=index,
        area=area,
        edge_lengths=(lab, lbc, lca),
        min_angle_deg=min(
            _angle_at(a, b, c), _angle_at(b, c, a), _angle_at(c, a, b)
        ),
        max_angle_deg=max(
            _angle_at(a, b, c), _angle_at(b, c, a), _angle_at(c, a, b)
        ),
        edge_ratio=long_ / short if short > 0.0 else math.inf,
        aspect_ratio=long_ * perimeter / (_FOUR_SQRT3 * area),
        aspect_frobenius=frobenius,
        gamma=1.0 / frobenius if frobenius > 0.0 else 0.0,
        radius_ratio=radius_ratio,
        scaled_jacobian=two_area / (mid * long_) if mid * long_ > 0.0 else 0.0,
        condition_number_inverse=sicn,
        degenerate=False,
    )


@dataclass(frozen=True)
class QualitySummary:
    """Aggregate of one measure over a mesh, with the element that produced the worst
    value -- the only element a repair tool needs to be pointed at."""

    metric: str
    minimum: float
    maximum: float
    mean: float
    standard_deviation: float
    worst_index: int
    worst_value: float
    comparison: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "metric": self.metric,
            "minimum": self.minimum,
            "maximum": self.maximum,
            "mean": self.mean,
            "standard_deviation": self.standard_deviation,
            "worst_index": self.worst_index,
            "worst_value": self.worst_value,
            "comparison": self.comparison,
        }


@dataclass(frozen=True)
class MeshQualityReport:
    """The whole measurement of a mesh: per-element records plus per-measure aggregates."""

    triangle_count: int
    degenerate_count: int
    summary: dict[str, QualitySummary]
    percentiles: dict[str, dict[str, float]] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return {
            "triangle_count": self.triangle_count,
            "degenerate_count": self.degenerate_count,
            "summary": {name: item.to_dict() for name, item in self.summary.items()},
            "percentiles": self.percentiles,
        }


def _percentile(sorted_values: Sequence[float], q: float) -> float:
    """Linear-interpolated percentile of an already sorted sequence (q in [0, 100])."""
    if not sorted_values:
        return 0.0
    if len(sorted_values) == 1:
        return sorted_values[0]
    position = (len(sorted_values) - 1) * q / 100.0
    low = math.floor(position)
    high = math.ceil(position)
    if low == high:
        return sorted_values[int(position)]
    weight = position - low
    return sorted_values[low] * (1.0 - weight) + sorted_values[high] * weight


def mesh_quality(
    mesh: Mesh,
    *,
    percentiles: Iterable[float] = (5.0, 50.0, 95.0),
    faces: Sequence[Sequence[int]] | None = None,
) -> MeshQualityReport:
    """Measure every face of ``mesh`` (or of ``faces``) and aggregate each measure.

    The returned report holds the min/max/mean/deviation of every measure *and* the
    index of the element that realises the gate-relevant extreme, so a caller can point
    a repair at ``summary["gamma"].worst_index`` instead of at the whole shell.
    """
    selected = mesh.faces if faces is None else faces
    records: list[TriangleQuality] = []
    for index, face in enumerate(selected):
        vertices = mesh.vertices
        records.append(
            triangle_quality(
                vertices[face[0]], vertices[face[1]], vertices[face[2]], index=index
            )
        )
    summary: dict[str, QualitySummary] = {}
    percentiles_out: dict[str, dict[str, float]] = {}
    for metric in METRIC_DEFINITIONS:
        comparison = METRIC_COMPARISON[metric]
        values = [record.value(metric) for record in records]
        finite = [value for value in values if math.isfinite(value)]
        if not values:
            summary[metric] = QualitySummary(metric, 0.0, 0.0, 0.0, 0.0, -1, 0.0, comparison)
            percentiles_out[metric] = {f"p{q:g}": 0.0 for q in percentiles}
            continue
        worst_index, worst_value = _extreme(values, comparison)
        mean = sum(finite) / len(finite) if finite else math.inf
        deviation = (
            math.sqrt(sum((value - mean) ** 2 for value in finite) / len(finite))
            if finite
            else math.inf
        )
        summary[metric] = QualitySummary(
            metric=metric,
            minimum=min(values) if finite else 0.0,
            maximum=max(values) if finite else math.inf,
            mean=mean,
            standard_deviation=deviation,
            worst_index=worst_index,
            worst_value=worst_value,
            comparison=comparison,
        )
        ordered = sorted(finite)
        percentiles_out[metric] = {f"p{q:g}": _percentile(ordered, q) for q in percentiles}
    return MeshQualityReport(
        triangle_count=len(records),
        degenerate_count=sum(1 for record in records if record.degenerate),
        summary=summary,
        percentiles=percentiles_out,
    )


def _extreme(values: Sequence[float], comparison: str) -> tuple[int, float]:
    """Index and value of the element the gate cares about (worst, not average)."""
    if comparison == "min":
        best_index = min(range(len(values)), key=lambda i: values[i])
    else:
        best_index = max(range(len(values)), key=lambda i: values[i])
    return best_index, values[best_index]


@dataclass(frozen=True)
class QualityCheck:
    """One gate decision, with the element that caused it and what to do about it."""

    metric: str
    comparison: str
    limit: float
    value: float
    passed: bool
    worst_index: int
    hint: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "metric": self.metric,
            "comparison": self.comparison,
            "limit": self.limit,
            "value": self.value,
            "passed": self.passed,
            "worst_index": self.worst_index,
            "hint": self.hint,
        }


@dataclass(frozen=True)
class QualityVerdict:
    """A CI door worth having: it says *which* measure failed, *where*, and *why*."""

    checks: tuple[QualityCheck, ...]
    report: MeshQualityReport

    @property
    def passed(self) -> bool:
        return all(check.passed for check in self.checks)

    @property
    def failures(self) -> tuple[QualityCheck, ...]:
        return tuple(check for check in self.checks if not check.passed)

    def exit_code(self) -> int:
        return 0 if self.passed else 1

    def to_dict(self) -> dict[str, Any]:
        return {
            "passed": self.passed,
            "exit_code": self.exit_code(),
            "checks": [check.to_dict() for check in self.checks],
            "report": self.report.to_dict(),
        }


def quality_gate(
    mesh: Mesh,
    *,
    criteria: Mapping[str, float] | None = None,
) -> QualityVerdict:
    """Run the quality gate on ``mesh``.

    ``criteria`` maps a measure of :data:`METRIC_DEFINITIONS` to its limit; the
    comparison sense comes from :data:`METRIC_COMPARISON`, never from the caller, so a
    limit cannot be inverted by hand.  ``None`` means :data:`DEFAULT_CRITERIA`.
    """
    active = dict(DEFAULT_CRITERIA if criteria is None else criteria)
    unknown = sorted(set(active) - set(METRIC_DEFINITIONS))
    if unknown:
        raise KeyError(f"unknown metric(s) in criteria: {', '.join(unknown)}")
    report = mesh_quality(mesh)
    checks: list[QualityCheck] = []
    for metric, limit in active.items():
        comparison = METRIC_COMPARISON[metric]
        item = report.summary[metric]
        value = item.minimum if comparison == "min" else item.maximum
        passed = value >= limit if comparison == "min" else value <= limit
        checks.append(
            QualityCheck(
                metric=metric,
                comparison=comparison,
                limit=float(limit),
                value=float(value),
                passed=passed,
                worst_index=item.worst_index,
                hint="" if passed else _REPAIR_HINTS.get(metric, ""),
            )
        )
    return QualityVerdict(checks=tuple(checks), report=report)


def gate_report(mesh: Mesh, *, criteria: Mapping[str, float] | None = None) -> dict[str, Any]:
    """JSON-safe convenience wrapper around :func:`quality_gate`."""
    return quality_gate(mesh, criteria=criteria).to_dict()


if __name__ == "__main__":  # pragma: no cover
    # Smoke path: an equilateral triangle at unit edge -- every ideal value at once.
    tri = triangle_quality((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.5, _SQRT3 / 2.0, 0.0))
    print(json.dumps(tri.to_dict(), indent=2))