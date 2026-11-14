"""CAD-side DG curved-wall hand-off: construction, measurement, refusal, frame.

The tests assert invariants (exactness against an analytic circle, digest
survival, refusal of malformed payloads), not merely that code runs.  The one
test that touches TensorDG is skipped when that package is not importable; it
is the cross-implementation check that the CAD evaluator and the consumer agree.
"""
from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path

import pytest

from tensorcad import dg_handoff as dh

ROOT = Path(__file__).resolve().parents[1]


def _ring_in_plane(radius: float, n: int = 32, offset_z: float = 0.0):
    return [(radius * math.cos(2 * math.pi * i / n),
             radius * math.sin(2 * math.pi * i / n), offset_z) for i in range(n)]


# --------------------------------------------------------------- construction


def test_circle_contract_is_the_exact_circle() -> None:
    """Tangent-intersection control points must reproduce a circle, not approximate it."""
    c = dh.circle_contract(1.0, (0.0, 0.0))
    r = dh.curve_geometry_report(c)
    assert r["max_radial_deviation_rel"] < 1e-14
    assert r["circumference_rel_error"] < 1e-14
    assert r["area_rel_error"] < 1e-14
    assert r["closure_gap_m"] == 0.0
    # rational quadratic circle: 2 segments + 1 control points per segment
    assert r["control_points"] == 9 and r["segments"] == 4 and r["degree"] == 2


def test_control_polygon_follows_the_cad_construction() -> None:
    """Middle control point = tangent intersection at ``r / cos(dtheta/2)``."""
    radius, segments = 2.0, 6
    c = dh.circle_contract(radius, (1.0, -0.5), segments=segments)
    cps, w, centre = c["control_points"], c["weights"], c["centre"]
    assert len(cps) == 2 * segments + 1 and len(w) == len(cps)
    assert cps[-1] == cps[0]
    cos_half = math.cos(math.pi / segments)
    for i in range(segments):
        corner = cps[2 * i + 1]
        d = math.hypot(corner[0] - centre[0], corner[1] - centre[1])
        assert d == pytest.approx(radius / cos_half, rel=1e-15)
        assert w[2 * i + 1] == pytest.approx(cos_half, rel=1e-15)
        on_circle = cps[2 * i]
        assert math.hypot(on_circle[0] - centre[0], on_circle[1] - centre[1]) == pytest.approx(
            radius, rel=1e-15)


def test_segment_count_is_a_representation_choice_not_a_geometry_choice() -> None:
    """A 4- and an 8-segment contract are the same circle in different splines."""
    r4 = dh.curve_geometry_report(dh.circle_contract(1.5, (0.3, -0.2), segments=4))
    r8 = dh.curve_geometry_report(dh.circle_contract(1.5, (0.3, -0.2), segments=8))
    assert r8["control_points"] == 17
    assert r8["circumference_m"] == pytest.approx(r4["circumference_m"], rel=1e-13)
    assert r8["enclosed_area_m2"] == pytest.approx(r4["enclosed_area_m2"], rel=1e-13)
    assert r8["max_radial_deviation_rel"] < 1e-14


def test_derivative_is_tangent_and_perpendicular_to_the_radius() -> None:
    c = dh.circle_contract(1.0)
    ts = [i / 97 for i in range(1, 97)]
    pts = dh.evaluate_curve(c, ts)
    der = dh.evaluate_curve_derivative(c, ts)
    for (x, y), (dx, dy) in zip(pts, der):
        assert abs(x * dx + y * dy) < 1e-13      # tangent . radius == 0
        assert math.hypot(dx, dy) > 1.0


def test_construction_rejects_nonsense() -> None:
    for kwargs in ({"radius": 0.0}, {"radius": -1.0}, {"radius": float("nan")}):
        with pytest.raises(dh.HandoffError):
            dh.circle_contract(**kwargs)
    with pytest.raises(dh.HandoffError):
        dh.circle_contract(1.0, segments=0)
    with pytest.raises(dh.HandoffError):
        dh.circle_contract(1.0, plane="yz")


# --------------------------------------------------------------------- frame


def test_plane_map_declares_chart_handedness() -> None:
    xy = dh.PlaneMap("xy", *dh.PLANES["xy"], 0.0)
    xz = dh.PlaneMap("xz", *dh.PLANES["xz"], -0.4)
    assert xy.chart_handedness == 1 and xz.chart_handedness == -1
    assert xz.project((1.0, 2.0, 3.0)) == (1.0, 3.0)      # (u, v) = (x, z)
    assert xz.out_of_plane((1.0, 2.0, 3.0)) == pytest.approx(2.0 + 0.4)
    assert xy.to_3d((5.0, 6.0)) == (5.0, 6.0, 0.0)
    assert dh.PlaneMap.from_dict(xz.to_dict()) == xz


def test_plane_check_catches_a_mirrored_chart() -> None:
    """The equator of a sphere may not be read as if it lived in the x-z plane."""
    xy_plane = dh.PlaneMap("xy", *dh.PLANES["xy"], 0.0)
    ring = _ring_in_plane(1.0)
    ok = dh.check_plane_ring(ring, xy_plane, radius=1.0, centre=(0.0, 0.0))
    assert ok["within_tolerance"] and ok["out_of_plane_max_m"] == 0.0
    bad = dh.check_plane_ring(ring, dh.PlaneMap("xz", *dh.PLANES["xz"], 0.0),
                              radius=1.0, centre=(0.0, 0.0))
    assert not bad["within_tolerance"]          # out-of-plane component ~1 m


# -------------------------------------------------------------------- bundle


def test_bundle_round_trip_and_digest(tmp_path: Path) -> None:
    ring = _ring_in_plane(1.0, 48)
    contract = dh.circle_contract(1.0)
    manifest = dh.write_dg_bundle(tmp_path / "w", contract, r_outer_m=2.0,
                                  source_3d_ring=ring, case="test")
    assert manifest["geometry"]["digest"] == dh.verify_bundle(tmp_path / "w")["digest"]
    check = dh.verify_bundle(tmp_path / "w")
    assert check["verified"] and check["r_outer_m"] == 2.0
    assert check["measured"]["max_radial_deviation_rel"] < 1e-14
    assert json.loads((tmp_path / "w" / "curved_wall.json").read_text())["degree"] == 2


def test_bundle_refuses_a_lying_handoff(tmp_path: Path) -> None:
    contract = dh.circle_contract(1.0)
    dh.write_dg_bundle(tmp_path / "w", contract, r_outer_m=2.0)
    handoff = tmp_path / "w" / "handoff.json"
    payload = json.loads(handoff.read_text())
    payload["measured"]["circumference_m"] *= 1.01
    handoff.write_text(json.dumps(payload))
    with pytest.raises(dh.HandoffError, match="circumference_m"):
        dh.verify_bundle(tmp_path / "w")


def test_bundle_refuses_a_tampered_contract(tmp_path: Path) -> None:
    contract = dh.circle_contract(1.0)
    dh.write_dg_bundle(tmp_path / "w", contract, r_outer_m=2.0)
    wall = tmp_path / "w" / "curved_wall.json"
    payload = json.loads(wall.read_text())
    payload["weights"][1] *= 1.02
    wall.write_text(json.dumps(payload))
    with pytest.raises(dh.HandoffError, match="digest"):
        dh.verify_bundle(tmp_path / "w")


def test_bundle_refuses_a_ring_off_the_declared_plane(tmp_path: Path) -> None:
    ring = [(x, y, 0.2) for x, y, _ in _ring_in_plane(1.0)]   # lifted off z = 0
    with pytest.raises(dh.HandoffError, match="does not project"):
        dh.write_dg_bundle(tmp_path / "w", dh.circle_contract(1.0), r_outer_m=2.0,
                           source_3d_ring=ring)


def test_bundle_refuses_an_outer_radius_inside_the_wall(tmp_path: Path) -> None:
    with pytest.raises(dh.HandoffError):
        dh.write_dg_bundle(tmp_path / "w", dh.circle_contract(1.0), r_outer_m=0.9)


def test_malformed_contract_is_refused() -> None:
    good = dh.circle_contract(1.0)
    for key, value in (("schema", "tensordg.curved_geometry.v0"), ("kind", "surface"),
                       ("degree", 0), ("closed", False)):
        bad = dict(good)
        bad[key] = value
        with pytest.raises(dh.HandoffError):
            dh.curve_geometry_report(bad)
    bad = dict(good)
    bad["knots"] = good["knots"][:-1]
    with pytest.raises(dh.HandoffError):
        dh.curve_geometry_report(bad)
    bad = dict(good)
    bad["weights"] = [0.0] + good["weights"][1:]
    with pytest.raises(dh.HandoffError):
        dh.curve_geometry_report(bad)


# ------------------------------------------------------- cross implementation


def test_matches_the_consumer_evaluator() -> None:
    """The CAD curve and TensorDG's de Boor evaluation must agree to round-off."""
    dg = pytest.importorskip("tensordg.curved_geometry")
    contract = dh.circle_contract(0.75, (-1.25, 2.0), segments=4)
    ts = [i / 400 for i in range(401)]
    mine = dh.evaluate_curve(contract, ts)
    theirs = dg.nurbs_curve_points(contract["degree"], contract["knots"],
                                   contract["control_points"], contract["weights"], ts)
    worst = max(math.hypot(a[0] - b[0], a[1] - b[1]) for a, b in zip(mine, theirs))
    assert worst < 1e-12
    dg.validate_geometry_contract(contract)     # the consumer's own gate accepts it


def test_module_is_standard_library_only() -> None:
    """The zero-dependency core lane must not need numpy/torch to build a contract."""
    import ast

    source = (ROOT / "tensorcad" / "dg_handoff.py").read_text()
    imported = set()
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, ast.Import):
            imported.update(alias.name.split(".")[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.add(node.module.split(".")[0])
    assert not imported & {"numpy", "torch", "scipy"}
