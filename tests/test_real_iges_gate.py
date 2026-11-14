"""Regression tests for the real-publisher IGES route (anchor GAP-CAD-T1-IGES-PIPE).

The defects these tests pin down were found by running the pipeline on the publisher's own
files (KCS / KVLCC2 / JBC / KP505 in the asset store), not on the self-written patch sets
the earlier rounds used:

* a seeded trim-loop inversion latched onto a domain wall and stayed pinned there, so the
  face was refused with an inversion residual of 10-17 file units;
* a *singular* edge -- a whole parametric edge whose 3-D image is one point -- inverted to
  a sliver and the ear clipper refused the face ("not a simple polygon");
* ear clipping only ever emits trim-loop vertices, so a patch had no interior sample at
  all (KCS wetted area came out 30 % low).
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from tensorcad.geomquality.trim_mesh import (
    _loop_to_uv,
    _rebuild_singular_edges,
    invert_point,
    refine_parametric,
    triangulate_trimmed_domain,
)

ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "assets" / "geometry_registry.json"


class CollapsedEdgeSurface:
    """u in [0,1], v in [0,1]; the entire ``u = 1`` edge maps to the single point (1,0,0)."""

    domain = ((0.0, 1.0), (0.0, 1.0))

    def evaluate(self, u: float, v: float) -> tuple[float, float, float]:
        return (u, (1.0 - u) * v, 0.0)


class CurvedPatch:
    """A quarter cylinder: constant curvature in u, straight in v."""

    domain = ((0.0, math.pi / 2.0), (0.0, 1.0))

    def evaluate(self, u: float, v: float) -> tuple[float, float, float]:
        return (math.cos(u), math.sin(u), v)


def _domain_rectangle(surface, samples: int = 8):
    (u0, u1), (v0, v1) = surface.domain
    points = [surface.evaluate(u0 + (u1 - u0) * i / samples, v0) for i in range(samples)]
    points += [surface.evaluate(u1, v0 + (v1 - v0) * i / samples) for i in range(samples)]
    points += [surface.evaluate(u1 - (u1 - u0) * i / samples, v1) for i in range(samples)]
    points += [surface.evaluate(u0, v1 - (v1 - v0) * i / samples) for i in range(samples)]
    return points


def test_inversion_finds_the_interior_parameter_of_a_curved_patch():
    """The Levenberg-damped step must not stall before it reaches an interior parameter."""
    surface = CurvedPatch()
    for u_true, v_true in ((0.3, 0.2), (1.1, 0.9), (0.0, 0.5)):
        point = surface.evaluate(u_true, v_true)
        u, v, residual = invert_point(surface, point, seeds=[(0.0, 0.0)],
                                      fallback_tolerance=1e-9)
        assert residual <= 1e-9, (u_true, residual)
        assert math.dist(surface.evaluate(u, v), point) <= 1e-9


def test_singular_edge_loop_stays_a_fillable_polygon():
    surface = CollapsedEdgeSurface()
    loop = _domain_rectangle(surface)
    uv, residual = _loop_to_uv(surface, loop, 1e-6)
    assert residual <= 1e-6
    triangles = triangulate_trimmed_domain([uv])
    assert len(triangles) >= 2
    area = sum(
        0.5 * abs((uv[b][0] - uv[a][0]) * (uv[c][1] - uv[a][1])
                  - (uv[b][1] - uv[a][1]) * (uv[c][0] - uv[a][0]))
        for a, b, c in triangles
    )
    assert abs(area - 1.0) < 0.05, area


def test_rebuild_singular_edges_spans_the_wall():
    wall = 1.0
    parameters = [(0.0, 0.0), (0.5, 0.0), (wall, 0.002), (wall, 0.002), (0.5, 1.0),
                  (0.0, 1.0)]
    points = [(0.0, 0.0, 0.0), (0.5, 0.0, 0.0), (wall, 0.0, 0.0), (wall, 0.0, 0.0),
              (0.5, 1.0, 0.0), (0.0, 1.0, 0.0)]
    rebuilt = _rebuild_singular_edges(parameters, points, (0.0, 1.0, 0.0, 1.0), 1e-9)
    assert (wall, 0.0) in rebuilt or (wall, pytest.approx(0.0)) in rebuilt
    assert any(abs(u - wall) < 1e-12 and abs(v - 1.0) < 1e-12 for u, v in rebuilt), rebuilt


def test_refinement_inserts_interior_samples():
    """Ear clipping emits trim-loop vertices only; refinement must add surface samples."""
    surface = CurvedPatch()
    uv = [(0.0, 0.0), (math.pi / 2.0, 0.0), (math.pi / 2.0, 1.0), (0.0, 1.0)]
    coarse = [(0, 1, 2), (0, 2, 3)]
    points, positions, triangles = refine_parametric(surface, uv, coarse, 0.02, levels=6)
    assert len(points) > len(uv), "no interior sample was inserted"
    assert len(triangles) > len(coarse)
    assert all(len({a, b, c}) == 3 for a, b, c in triangles), "degenerate triangle"
    # boundary vertices keep their exact parameters, so shared seams are unchanged
    assert points[:4] == uv
    # every new sample lies on the surface by construction
    for index in range(len(uv), len(points)):
        u, v = points[index]
        assert math.dist(positions[index], surface.evaluate(u, v)) <= 1e-12


def _asset(name: str) -> Path | None:
    if not REGISTRY.exists():
        return None
    registry = json.loads(REGISTRY.read_text())
    store = Path(registry.get("store_root", ""))
    for asset in registry["assets"]:
        if asset["asset_id"].endswith(name):
            path = store / asset["asset_id"][len("asset://"):] if store else None
            return path if path and path.exists() else None
    return None


def test_publisher_kcs_hull_imports_without_a_refused_face():
    """The real KCS specimen: every patch must tessellate and the beam must be right."""
    source = _asset("KCS_hull_SVA.igs")
    if source is None:
        pytest.skip("publisher geometry not present in the asset store")
    from tensorcad.geomquality.trim_mesh import import_cad

    result = import_cad(source)
    report = result.report
    assert report["import"]["faces_in_file"] == 55
    assert report["import"]["faces_refused"] == 0, report["import"]["refused"]
    assert report["import"]["max_inversion_residual"] <= report["import"]["residual_tolerance"]
    assert report["before_repair"]["triangles"] > 3000
    vertices = result.mesh.vertices
    beam = 2.0 * max(point[1] for point in vertices)
    assert abs(beam - 32.2) / 32.2 < 0.01, beam
