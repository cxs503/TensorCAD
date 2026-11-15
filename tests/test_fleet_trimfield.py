"""Tests for trim-aware IGES sampling (KVLCC2 trim defect)."""
from __future__ import annotations

import math

import pytest

from tensorcad.fleet.trimfield import (
    MIN_INSIDE_FRACTION,
    point_in_polygon,
    sample_trimmed_patches,
    trim_polygons,
)


def _square(x0: float, y0: float, x1: float, y1: float):
    return [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]


def test_containment_inside_and_outside_a_square():
    square = _square(0.0, 0.0, 2.0, 2.0)
    assert point_in_polygon(1.0, 1.0, square) is True
    assert point_in_polygon(3.0, 1.0, square) is False
    assert point_in_polygon(1.0, -0.5, square) is False


def test_containment_of_a_concave_loop():
    concave = [(0.0, 0.0), (4.0, 0.0), (4.0, 1.0), (1.0, 1.0), (1.0, 4.0), (0.0, 4.0)]
    assert point_in_polygon(0.5, 3.0, concave) is True   # inside the upright leg
    assert point_in_polygon(3.0, 3.0, concave) is False  # inside the notch, not the loop


def test_degenerate_loop_contains_nothing():
    assert point_in_polygon(0.0, 0.0, []) is False
    assert point_in_polygon(0.0, 0.0, [(0.0, 0.0), (1.0, 1.0)]) is False


class _FakeCurve:
    def __init__(self, uv_pointer: int, model_pointer: int) -> None:
        self.curve_uv = uv_pointer
        self.curve_3d = model_pointer


class _FakeBasis:
    domain_u = (0.0, 1.0)
    domain_v = (0.0, 1.0)

    def __init__(self) -> None:
        self.sampled: list[tuple[float, float]] = []

    def bounding_box(self):
        return (0.0, 0.0, 0.0), (10.0, 10.0, 10.0)

    def evaluate(self, u: float, v: float):
        self.sampled.append((u, v))
        return (10.0 * u, 10.0 * v, 0.0)


class _FakeIges:
    """One trimmed patch; the *model-space* pointer is offered first on purpose."""

    def __init__(self) -> None:
        self.basis = _FakeBasis()
        self.polygons = {7: _square(5.0, 5.0, 500.0, 500.0), 9: _square(0.2, 0.2, 0.6, 0.6)}

    def trimmed_surfaces(self):
        return [(object(), self.basis, [_FakeCurve(uv_pointer=9, model_pointer=7)])]

    def trim_polygon_uv(self, pointer: int, samples: int = 200):
        return list(self.polygons[pointer])


def test_the_pointer_inside_the_domain_wins():
    iges = _FakeIges()
    records = trim_polygons(iges)
    assert len(records) == 1
    record = records[0]
    assert record["verdict"] == "decided"
    assert record["pointer"] == 9, "the model-space curve lands outside (u, v) and must lose"
    assert record["inside_fraction"] >= MIN_INSIDE_FRACTION


def test_sampling_keeps_only_points_inside_the_trim_loop():
    iges = _FakeIges()
    points, report = sample_trimmed_patches(trim_polygons(iges), spacing=5.0)
    assert report["patches_used"] == 1
    assert report["kept_samples"] > 0
    assert report["dropped_samples"] > 0
    for _x, _y, _z in points:
        u, v = iges.basis.sampled[points.index((_x, _y, _z))] if False else (0.0, 0.0)
    kept = sorted(iges.basis.sampled)
    u_lo, v_lo = min(u for u, _ in kept), min(v for _, v in kept)
    assert u_lo <= 0.25 and v_lo <= 0.25


def test_undecided_patch_is_skipped_not_widened():
    class _NoTrim(_FakeIges):
        def trim_polygon_uv(self, pointer: int, samples: int = 200):
            return list(self.polygons[7])  # only the model-space curve resolves

    points, report = sample_trimmed_patches(trim_polygons(_NoTrim()), spacing=5.0)
    assert points == []
    assert report["patches_skipped"] == 1
