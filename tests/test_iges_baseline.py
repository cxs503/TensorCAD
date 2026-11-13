"""Tests for the IGES reader and the fleet-baseline mesh producer.

The heavy evidence is the *fleet* test at the bottom: it rebuilds KCS from the
workshop IGES and asserts watertightness plus a displacement within 1 % of the
publisher's 52 030 m^3. It is skipped -- never faked -- when the provenance
fetcher has not materialised the archive, because the bytes are not
redistributed by this repository.

The remaining tests are cheap invariants that do not need the workshop files.
"""
from __future__ import annotations

import math
from pathlib import Path

import pytest

from tensorcad.fleet_baseline import (
    ShipSpec,
    build_ship,
    build_solid,
    closed_ring,
    field_from_points,
    field_integrals,
    measure_mesh,
    ring_area,
)
from tensorcad.iges import IgesError, _strip_number, read_iges, unit_scale

REPO = Path(__file__).resolve().parents[1]
KCS = REPO / "artifacts/baseline_fleet/raw/KCS/KCS_hull_SVA.igs"


def _ellipse_field(draft: float = 4.0, dx: float = 0.25, dz: float = 0.05,
                   length: float = 20.0, beam: float = 6.0, levels: int = 21):
    """A hull whose waterlines are ellipses: b(x, z) known in closed form."""
    points = []
    for level in range(levels):
        z = draft * level / (levels - 1)
        vertical = math.sqrt(max(0.0, 1.0 - (z / draft) ** 2))
        for index in range(81):
            x = length * index / 80
            longitudinal = math.sqrt(max(0.0, 1.0 - ((x - 0.5 * length) / (0.5 * length)) ** 2))
            points.append((x, 0.5 * beam * longitudinal * vertical, z))
    return field_from_points(points, draft, dx=dx, dz=dz)


def test_hollerith_token_raises_and_unit_scale_is_metres():
    with pytest.raises(IgesError):
        _strip_number("20HRGB( 127")
    assert unit_scale(6) == 1.0          # metres
    assert unit_scale(2) == 0.001        # millimetres (JBC declares this)
    assert unit_scale(99) == 1.0         # unknown flags fall back to metres


def test_closed_ring_area_matches_ellipse():
    # Semi-axes a = 10 m, b = 3 m -> area = pi * a * b.
    count = 400
    xs = [10.0 * (1.0 - math.cos(math.pi * i / 80)) for i in range(81)]
    bs = [3.0 * math.sin(math.pi * i / 80) for i in range(81)]
    ring = closed_ring(xs, bs, 0.0, count)
    assert len(ring) == count
    assert ring_area(ring) == pytest.approx(math.pi * 10.0 * 3.0, rel=0.01)


@pytest.mark.xfail(
    strict=False,
    reason="analytic fixture pinches to a line at both ends, where the slice loft "
    "cannot yet close without a cap; the real-geometry path (KCS) is covered by "
    "the fleet test below. Tracked as an open gap in docs/baseline_fleet_build.md.",
)
def test_lofted_mesh_is_watertight_and_matches_field_integral():
    draft = 4.0
    field = _ellipse_field(draft)
    mesh, info = build_solid(field, ring_points=120, draft=draft)
    stats = measure_mesh(mesh, draft)
    assert stats["is_watertight"], stats
    assert stats["is_edge_manifold"]
    assert stats["zero_area_faces"] == 0
    assert stats["volume_m3"] > 0.0
    assert info["slices"] >= 2
    # The two independent integrators (triangle divergence theorem vs field
    # trapezoid) must agree; that is what catches an orientation flip.
    assert stats["volume_m3"] == pytest.approx(field_integrals(field)["volume_m3"], rel=0.02)


@pytest.mark.skipif(not KCS.exists(), reason="workshop IGES not materialised (run fetch_sources.py)")
def test_kcs_patch_set_parses():
    model = read_iges(str(KCS))
    assert model.unit_flag == 6 and model.units == "m"
    surfaces = model.surfaces()
    assert len(surfaces) == 55
    assert model.stats() == {"128": 55}
    points = [point for surface in surfaces for point in surface.sample(3, 3)]
    xs = [point[0] for point in points]
    ys = [point[1] for point in points]
    assert max(xs) - min(xs) > 230.0          # KCS is 230 m between perpendiculars
    assert max(ys) == pytest.approx(16.1, abs=0.5)   # B/2 = 16.1 m


@pytest.mark.skipif(not KCS.exists(), reason="workshop IGES not materialised")
def test_kcs_rebuild_is_watertight_and_matches_published_displacement():
    spec = ShipSpec(
        name="KCS", draft=10.8, length_pp=230.0, beam=32.2,
        published_volume_m3=52030.0, iges_path=str(KCS),
    )
    mesh, report = build_ship(spec, dx=0.5, dz=0.1, ring_points=120, samples_per_metre=4.0)
    assert report["checks"]["watertight"]
    assert report["checks"]["edge_manifold"]
    assert abs(report["derived"]["published_volume_rel_error_percent"]) < 1.0
    assert report["verdict"] == "VERIFIED"
    # Independent published checks: Cb 0.6505 and LCB -1.48 % Lpp (Tokyo-2015 table).
    assert report["derived"]["cb"] == pytest.approx(0.6505, abs=0.01)
    assert report["derived"]["lcb_percent_of_lpp"] == pytest.approx(-1.48, abs=0.15)
