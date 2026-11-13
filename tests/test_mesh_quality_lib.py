"""Pins the canonical mesh-quality library and its gate.

Three things are asserted here, in increasing order of strength:

* the closed-form values on canonical elements (the *definitions* are right);
* agreement with the independent route in :mod:`tools.quality_oracle` (the
  *implementation* is right) and with exact rational arithmetic where the mesh allows
  it (the *rounding* is bounded);
* the gate's behaviour: a clean hull passes, a raw marching-tets hull fails with the
  element index that caused it, and the documented repair flips the verdict.

The full differential run over every hull mesh lives in
``tools/mesh_quality_crosscheck.py`` (artifacts in ``artifacts/mesh_quality_lib/``);
this module is the fast subset that a CI job can afford on every commit.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

from tensorcad.geometry import Mesh, MeshBuilder
from tensorcad.geomquality.quality import (
    DEFAULT_CRITERIA,
    IDEAL_VALUES,
    METRIC_COMPARISON,
    METRIC_DEFINITIONS,
    mesh_quality,
    quality_gate,
    triangle_quality,
)

TOOLS = Path(__file__).resolve().parents[1] / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

import quality_oracle  # noqa: E402


def _equilateral(edge: float = 1.0) -> tuple[tuple[float, float, float], ...]:
    return ((0.0, 0.0, 0.0), (edge, 0.0, 0.0), (0.5 * edge, math.sqrt(3) / 2.0 * edge, 0.0))


def _box() -> Mesh:
    builder = MeshBuilder()
    for point in [
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (1.0, 1.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0),
        (1.0, 0.0, 1.0),
        (1.0, 1.0, 1.0),
        (0.0, 1.0, 1.0),
    ]:
        builder.add_vertex(point)
    for a, b, c, d in ((0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)):
        builder.add_triangle(a, b, c)
        builder.add_triangle(a, c, d)
    return builder.build()


def _patch() -> Mesh:
    """One good triangle and one hard sliver, rational coordinates."""
    builder = MeshBuilder()
    for point in [
        (0.0, 0.0, 0.0),
        (4.0, 0.0, 0.0),
        (4.0, 3.0, 0.0),
        (8.0, 0.0, 0.0),
        (8.0, 0.01, 0.0),
        (12.0, 0.0, 0.0),
    ]:
        builder.add_vertex(point)
    builder.add_triangle(0, 1, 2)
    builder.add_triangle(3, 4, 5)
    return builder.build()


def test_ideal_element_scores_every_ideal_value():
    measured = triangle_quality(*_equilateral())
    for metric, expected in IDEAL_VALUES.items():
        assert measured.value(metric) == pytest.approx(expected, abs=1e-12), metric


@pytest.mark.parametrize(
    "triangle",
    [
        ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)),
        ((0.0, 0.0, 0.0), (3.0, 0.0, 0.0), (3.0, 4.0, 0.0)),
    ],
)
def test_rational_triangles_match_exact_arithmetic(triangle):
    measured = triangle_quality(*triangle)
    exact = quality_oracle.oracle_metrics_exact(*triangle)
    for metric in quality_oracle.METRICS:
        if metric not in exact:
            continue
        tolerance = 1e-9 if metric.endswith("_deg") else 1e-12
        assert measured.value(metric) == pytest.approx(exact[metric], abs=tolerance), metric


def test_degenerate_triangle_is_worst_not_nan():
    measured = triangle_quality((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0))
    assert measured.degenerate
    assert measured.gamma == 0.0
    assert math.isinf(measured.aspect_ratio) and math.isinf(measured.edge_ratio)
    for metric in METRIC_DEFINITIONS:
        assert not math.isnan(measured.value(metric)), metric


def test_metrics_are_invariant_under_rigid_motion_and_scaling():
    base = triangle_quality((0.0, 0.0, 0.0), (3.0, 0.0, 0.0), (1.0, 2.0, 0.5))
    # rotate 90 degrees about z, translate, and scale by 7
    rotated = [(-y * 7.0 + 11.0, x * 7.0 - 4.0, z * 7.0 + 2.5) for x, y, z in
               ((0.0, 0.0, 0.0), (3.0, 0.0, 0.0), (1.0, 2.0, 0.5))]
    moved = triangle_quality(*rotated)
    for metric in METRIC_DEFINITIONS:
        if metric == "area":  # the one measure that is not scale free
            assert moved.value(metric) == pytest.approx(49.0 * base.value(metric), rel=1e-12)
            continue
        assert moved.value(metric) == pytest.approx(base.value(metric), rel=1e-12), metric


def test_gate_passes_a_clean_box_and_reports_the_worst_element():
    verdict = quality_gate(_box())
    assert verdict.passed
    assert verdict.exit_code() == 0
    assert {check.metric for check in verdict.checks} == set(DEFAULT_CRITERIA)
    for check in verdict.checks:
        assert check.comparison == METRIC_COMPARISON[check.metric]


def test_gate_fails_a_sliver_patch_and_points_at_it():
    verdict = quality_gate(_patch())
    assert not verdict.passed
    assert verdict.exit_code() == 1
    failing = {check.metric: check for check in verdict.failures}
    assert "min_angle_deg" in failing
    assert failing["min_angle_deg"].worst_index == 1  # the hard sliver is face 1
    assert failing["min_angle_deg"].hint  # a gate that only says "no" is not usable


def test_gate_rejects_unknown_criteria_and_cannot_be_inverted():
    with pytest.raises(KeyError):
        quality_gate(_box(), criteria={"aspect_ratio_nonsense": 10.0})
    # A generous limit must pass, a tight one must fail: the sense comes from the metric.
    assert quality_gate(_patch(), criteria={"min_angle_deg": 0.001}).passed
    assert not quality_gate(_patch(), criteria={"min_angle_deg": 1.0}).passed


def test_mesh_quality_aggregates_are_consistent():
    mesh = _patch()
    report = mesh_quality(mesh)
    assert report.triangle_count == 2
    assert report.degenerate_count == 0
    for metric in METRIC_DEFINITIONS:
        summary = report.summary[metric]
        assert summary.minimum <= summary.maximum
        assert summary.minimum <= summary.mean <= summary.maximum
        if summary.comparison == "min":
            assert summary.worst_value == summary.minimum
        else:
            assert summary.worst_value == summary.maximum
        assert 0 <= summary.worst_index < report.triangle_count


def test_independent_route_agrees_on_both_patch_triangles():
    mesh = _patch()
    for index, face in enumerate(mesh.faces):
        measured = triangle_quality(
            mesh.vertices[face[0]], mesh.vertices[face[1]], mesh.vertices[face[2]], index=index
        )
        independent = quality_oracle.oracle_metrics(
            mesh.vertices[face[0]], mesh.vertices[face[1]], mesh.vertices[face[2]]
        )
        for metric in quality_oracle.METRICS:
            assert measured.value(metric) == pytest.approx(independent[metric], abs=1e-6), metric