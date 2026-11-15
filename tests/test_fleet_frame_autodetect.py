"""Test for coordinate-frame autodetection in field_from_points (m21409).

The KCS/KVLCC2/JBC/ONRT official IGES files use a **keel-at-origin** frame
(z ∈ [0, draft], waterline at z=draft), while the legacy surrogate geometry
uses a **waterline-at-origin** frame (z ∈ [-draft, 0]).  The fix in m21409
makes ``field_from_points`` detect the frame automatically from the point
cloud's z-extent.

This test verifies both frames produce correct fields.
"""
from __future__ import annotations

import math

import pytest

from tensorcad.fleet_baseline import field_from_points


def _keel_at_origin_points(draft=10.0, beam=20.0, length=100.0, n=20):
    """Synthetic hull points in keel-at-origin frame (z ∈ [0, draft])."""
    pts = []
    for i in range(n):
        xi = i / (n - 1)
        x = length * xi
        taper = 1.0 - (2 * xi - 1) ** 2
        for k in range(n):
            z = draft * k / (n - 1)
            b = 0.5 * beam * taper * math.sqrt(z / draft) if z > 0 else 0.0
            if b > 1e-9:
                pts.append((x, b, z))
                pts.append((x, -b, z))
    return pts


def _waterline_at_origin_points(draft=10.0, beam=20.0, length=100.0, n=20):
    """Synthetic hull points in waterline-at-origin frame (z ∈ [-draft, 0])."""
    pts = []
    for i in range(n):
        xi = i / (n - 1)
        x = length * xi
        taper = 1.0 - (2 * xi - 1) ** 2
        for k in range(n):
            z = -draft * (1 - k / (n - 1))
            b = 0.5 * beam * taper * math.sqrt(-z / draft) if z < 0 else 0.0
            if b > 1e-9:
                pts.append((x, b, z))
                pts.append((x, -b, z))
    return pts


def test_keel_at_origin_frame():
    """Keel-at-origin (z ∈ [0, draft]) should produce a field with data."""
    draft = 10.0
    pts = _keel_at_origin_points(draft=draft)
    field = field_from_points(pts, draft=draft, dx=5.0, dz=0.5)
    assert field.nz > 0
    assert field.filled_nodes > 0
    # Max z with data should be at or near the waterline
    max_z = field.to_dict()["max_z_with_data_m"]
    assert abs(max_z - draft) < 1.0, f"max_z={max_z} should be near draft={draft}"


def test_waterline_at_origin_frame():
    """Waterline-at-origin (z ∈ [-draft, 0]) should also produce a correct field."""
    draft = 10.0
    pts = _waterline_at_origin_points(draft=draft)
    field = field_from_points(pts, draft=draft, dx=5.0, dz=0.5)
    assert field.nz > 0
    assert field.filled_nodes > 0
    # Should have most levels with data (hull tapers at ends)
    assert field.to_dict()["levels_with_data"] >= field.nz - 2


def test_keel_at_origin_max_breadth():
    """Max half-breadth should match the input beam/2."""
    draft = 10.0
    beam = 20.0
    pts = _keel_at_origin_points(draft=draft, beam=beam)
    field = field_from_points(pts, draft=draft, dx=5.0, dz=0.5)
    max_b = field.to_dict()["max_half_breadth_m"]
    assert abs(max_b - beam / 2) < 1.0


def test_no_points_raises():
    """Empty point cloud should raise ValueError."""
    with pytest.raises(ValueError, match="no surface samples"):
        field_from_points([], draft=10.0)


def test_points_above_waterline_only_raises():
    """Points entirely above the waterline should raise ValueError."""
    pts = [(50.0, 10.0, 15.0), (50.0, -10.0, 15.0)]
    with pytest.raises(ValueError, match="no surface samples in the immersed"):
        field_from_points(pts, draft=10.0)
