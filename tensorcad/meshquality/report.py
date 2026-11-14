"""Aggregate statistics and a machine gate on top of :mod:`tensorcad.meshquality.metrics`.

A metric library that only prints a minimum is not usable in CI: a pipeline needs to
say *which* threshold, on *which* cell kind, failed by *how much*, and it needs that
answer as data. This module is that layer:

* :func:`measure` turns points plus connectivity into a :class:`MeshQuality`;
* :meth:`MeshQuality.summary` reduces per-element values to statistics that survive
  outliers (minimum, 1st/5th percentile, median, mean, maximum) plus the count of
  non-positive (inverted or flat) cells;
* :class:`QualityGate` turns thresholds into :class:`Check` records with pass/fail,
  and :meth:`QualityGate.evaluate` returns a :class:`GateReport` that serialises to
  JSON. The CLI exits non-zero on a failed gate, which is the CI contract.

Default thresholds are *policy*, not measurements: they are stated in
:data:`DEFAULT_LIMITS`, overridable per call, and the benchmark example prints the
measured margin on every mesh so a reader can see how much room the policy leaves.
"""
from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from .metrics import (
    ELEMENT_KINDS,
    NODES_PER_CELL,
    ElementQuality,
    hex_metrics,
    prism_metrics,
    pyramid_metrics,
    quad_metrics,
    tet_metrics,
    tri_metrics,
)

__all__ = [
    "DEFAULT_GATE",
    "Check",
    "GateReport",
    "MeshQuality",
    "MetricLimit",
    "QualityGate",
    "SummaryStats",
    "measure",
    "summary_stats",
]

Vec3 = tuple[float, float, float]

_MEASURERS = {
    "tri": tri_metrics,
    "quad": quad_metrics,
    "tet": tet_metrics,
    "hex": hex_metrics,
    "prism": prism_metrics,
    "pyramid": pyramid_metrics,
}


# --------------------------------------------------------------------------------------
# statistics
# --------------------------------------------------------------------------------------
@dataclass(frozen=True)
class SummaryStats:
    """Outlier-resistant statistics of one metric over one cell kind."""

    metric: str
    count: int
    minimum: float
    p01: float
    p05: float
    median: float
    mean: float
    maximum: float
    non_positive: int

    def to_dict(self) -> dict[str, Any]:
        return {
            "metric": self.metric,
            "count": self.count,
            "min": self.minimum,
            "p01": self.p01,
            "p05": self.p05,
            "median": self.median,
            "mean": self.mean,
            "max": self.maximum,
            "non_positive": self.non_positive,
        }


def _percentile(sorted_values: Sequence[float], q: float) -> float:
    """Linear-interpolation percentile, same definition as numpy's default."""
    if not sorted_values:
        return math.nan
    if len(sorted_values) == 1:
        return sorted_values[0]
    position = q * (len(sorted_values) - 1)
    lower = math.floor(position)
    upper = min(lower + 1, len(sorted_values) - 1)
    weight = position - lower
    return sorted_values[lower] * (1.0 - weight) + sorted_values[upper] * weight


def summary_stats(metric: str, values: Iterable[float]) -> SummaryStats:
    data = [v for v in values if not math.isnan(v)]
    if not data:
        nan = math.nan
        return SummaryStats(metric, 0, nan, nan, nan, nan, nan, nan, 0)
    ordered = sorted(data)
    return SummaryStats(
        metric=metric,
        count=len(ordered),
        minimum=ordered[0],
        p01=_percentile(ordered, 0.01),
        p05=_percentile(ordered, 0.05),
        median=_percentile(ordered, 0.5),
        mean=sum(ordered) / len(ordered),
        maximum=ordered[-1],
        non_positive=sum(1 for v in ordered if v <= 0.0),
    )


# --------------------------------------------------------------------------------------
# container
# --------------------------------------------------------------------------------------
@dataclass(frozen=True)
class MeshQuality:
    """Per-element metrics of every cell kind present in one mesh."""

    elements: dict[str, ElementQuality] = field(default_factory=dict)

    def __getitem__(self, kind: str) -> ElementQuality:
        return self.elements[kind]

    def summary(self) -> dict[str, dict[str, SummaryStats]]:
        return {
            kind: {name: summary_stats(name, values) for name, values in quality.metrics.items()}
            for kind, quality in self.elements.items()
        }

    def to_dict(self) -> dict[str, Any]:
        return {
            kind: {
                "count": len(quality),
                "summary": {name: stat.to_dict() for name, stat in self.summary()[kind].items()},
            }
            for kind, quality in self.elements.items()
        }


def measure(
    points: Sequence[Sequence[float]],
    elements: Mapping[str, Sequence[Sequence[int]]],
) -> MeshQuality:
    """Measure every cell kind in ``elements``; unknown kinds are a ``ValueError``."""
    vertices: list[Vec3] = [(float(p[0]), float(p[1]), float(p[2])) for p in points]
    measured: dict[str, ElementQuality] = {}
    for kind, cells in elements.items():
        if kind not in _MEASURERS:
            raise ValueError(f"unknown cell kind {kind!r}; known kinds: {', '.join(ELEMENT_KINDS)}")
        if not cells:
            continue
        worst = max(len(cell) for cell in cells)
        expected = NODES_PER_CELL[kind]
        if worst != expected:
            raise ValueError(f"{kind} cells need {expected} nodes, got a cell with {worst}")
        for cell in cells:
            if any(i < 0 or i >= len(vertices) for i in cell):
                raise ValueError(f"{kind} cell {tuple(cell)} indexes outside {len(vertices)} points")
        measured[kind] = _MEASURERS[kind](vertices, cells)
    if not measured:
        raise ValueError("no cells to measure")
    return MeshQuality(measured)


# --------------------------------------------------------------------------------------
# gate
# --------------------------------------------------------------------------------------
@dataclass(frozen=True)
class Check:
    """One threshold applied to one metric of one cell kind."""

    name: str
    kind: str
    metric: str
    passed: bool
    value: float
    limit: str
    detail: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "kind": self.kind,
            "metric": self.metric,
            "passed": self.passed,
            "value": self.value,
            "limit": self.limit,
            "detail": self.detail,
        }


@dataclass(frozen=True)
class MetricLimit:
    """A lower and/or upper bound on one metric, plus how many bad cells are tolerated."""

    metric: str
    minimum: float | None = None
    maximum: float | None = None
    max_non_positive: int = 0
    note: str = ""

    def check(self, kind: str, stats: SummaryStats) -> Check:
        problems: list[str] = []
        if self.minimum is not None and not (stats.minimum >= self.minimum):
            problems.append(f"min {stats.minimum:.6g} < {self.minimum:g}")
        if self.maximum is not None and not (stats.maximum <= self.maximum):
            problems.append(f"max {stats.maximum:.6g} > {self.maximum:g}")
        if stats.non_positive > self.max_non_positive:
            problems.append(f"{stats.non_positive} non-positive cells > {self.max_non_positive}")
        limit = []
        if self.minimum is not None:
            limit.append(f">= {self.minimum:g}")
        if self.maximum is not None:
            limit.append(f"<= {self.maximum:g}")
        limit.append(f"non-positive <= {self.max_non_positive}")
        return Check(
            name=f"{kind}.{self.metric}",
            kind=kind,
            metric=self.metric,
            passed=not problems,
            value=stats.minimum if self.minimum is not None else stats.maximum,
            limit=" ".join(limit),
            detail=(
                "; ".join(problems)
                if problems
                else f"min {stats.minimum:.6g}, max {stats.maximum:.6g} over {stats.count} cells"
            ),
        )


@dataclass(frozen=True)
class QualityGate:
    """Thresholds per cell kind. A metric with no limit is reported, never gated."""

    limits: Mapping[str, tuple[MetricLimit, ...]]

    def evaluate(self, quality: MeshQuality) -> "GateReport":
        summary = quality.summary()
        checks: list[Check] = []
        for kind, limits in self.limits.items():
            if kind not in summary:
                continue  # a gate on a cell kind the mesh does not contain is vacuous
            for limit in limits:
                stats = summary[kind].get(limit.metric)
                if stats is None:
                    raise ValueError(f"gate names metric {limit.metric!r}, missing from {kind}")
                checks.append(limit.check(kind, stats))
        return GateReport(tuple(checks))


@dataclass(frozen=True)
class GateReport:
    checks: tuple[Check, ...]

    @property
    def passed(self) -> bool:
        return all(check.passed for check in self.checks)

    @property
    def failures(self) -> tuple[Check, ...]:
        return tuple(check for check in self.checks if not check.passed)

    def to_dict(self) -> dict[str, Any]:
        return {
            "passed": self.passed,
            "checks": [check.to_dict() for check in self.checks],
            "failures": [check.name for check in self.failures],
        }

    def write(self, path: str | Path) -> Path:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(self.to_dict(), indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return target


#: Default gate. The numbers are policy (documented in ``docs/mesh_quality.md``), not
#: measured thresholds: ``min_sicn`` 0.1 is roughly Gmsh's "the optimiser has work to do"
#: band, ``min_angle`` 10 deg is the usual surface-tessellation floor, and
#: ``edge_ratio`` 50 keeps a boundary-layer cell from being an accident. Every limit is
#: overridable; the benchmark reports the measured margin against them.
DEFAULT_LIMITS: dict[str, tuple[MetricLimit, ...]] = {
    "tri": (
        MetricLimit("min_sicn", minimum=0.1, note="shape and anisotropy"),
        MetricLimit("scaled_jacobian", minimum=0.1, note="shape only"),
        MetricLimit("min_angle", minimum=10.0, note="degrees"),
        MetricLimit("edge_ratio", maximum=50.0),
    ),
    "quad": (
        MetricLimit("min_sicn", minimum=0.1),
        MetricLimit("scaled_jacobian", minimum=0.1),
        MetricLimit("min_angle", minimum=10.0),
    ),
    "tet": (
        MetricLimit("min_sicn", minimum=0.1),
        MetricLimit("scaled_jacobian", minimum=0.0),
        MetricLimit("min_dihedral", minimum=5.0, note="degrees"),
    ),
    "hex": (
        MetricLimit("min_sicn", minimum=0.1),
        MetricLimit("scaled_jacobian", minimum=0.1),
        MetricLimit("min_angle", minimum=10.0),
    ),
    # A boundary-layer prism is *meant* to be anisotropic, so ``edge_ratio`` is not gated
    # (height vs base edge is the point of the layer); ``min_sicn`` catches the distortion
    # that actually hurts the discretisation, and it is aspect sensitive.
    "prism": (
        MetricLimit("min_sicn", minimum=0.1, note="shape and anisotropy"),
        MetricLimit("scaled_jacobian", minimum=0.1, note="shape only; ideal prism = 0.866"),
        MetricLimit("min_angle", minimum=10.0, note="degrees"),
    ),
    # Ideal-pyramid skewness is 0.1968, not 0 (a square pyramid has no 60 deg apex angle),
    # so skewness is not gated here; see docs/mesh_quality_prism_pyramid.md.
    "pyramid": (
        MetricLimit("min_sicn", minimum=0.1),
        MetricLimit("scaled_jacobian", minimum=0.0),
        MetricLimit("min_angle", minimum=10.0),
    ),
}

DEFAULT_GATE = QualityGate(DEFAULT_LIMITS)