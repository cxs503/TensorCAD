"""Pins the ADIABATIC relation between the two aspect-ratio dialects in this repo.

Audit context (m17437, anchor GAP-CAD-D1-MESH-QUALITY-LIB):

``tensorcad.mesh_conditioning`` (and ``tensorcad.validation``) have shipped
``aspect = l_max^2 / (2A)`` for a long time, with a production limit of 25.
``tensorcad.geomquality.quality`` ships ``aspect = l_max*p / (4*sqrt(3)*A)``,
also with a default limit of 25, and its docstring claims the two limits are
"the same measure re-expressed".

They are not the same measure: the ratio is ``p / (2*sqrt(3)*l_max)``, which
lies in ``(1/sqrt(3), sqrt(3)/2]`` -- always below 1, but not a constant. So a
threshold of 25 does not transfer between the dialects, and an element that the
old gate rejects can pass the new one.

These tests do not "fix" the dialects; they make the conversion explicit and
machine-checked so that no future edit can silently redefine either side and
leave the docstring claim looking true.
"""
from __future__ import annotations

import math

import pytest

# The audited module ships in m17404's branch; until it is merged into this branch the
# tests skip instead of failing, and they activate automatically the day it lands.
pytest.importorskip("tensorcad.geomquality.quality")
from tensorcad.geomquality.quality import triangle_quality  # noqa: E402

#: Two dialects of "aspect ratio" that must never be conflated.
LEGACY_LIMIT = 25.0
LIBRARY_LIMIT = 25.0

TRIS = {
    "equilateral": ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.5, math.sqrt(3.0) / 2.0, 0.0)),
    "right_3_4_5": ((0.0, 0.0, 0.0), (4.0, 0.0, 0.0), (0.0, 3.0, 0.0)),
    "sliver_2_2_0.02": ((0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (1.0, 0.02, 0.0)),
}


def legacy_aspect(tri) -> float:
    """``mesh_conditioning._aspect_ratio``: l_max^2 / |u x v| (we keep the helper local
    so the test does not couple to a private name)."""
    a, b, c = tri
    longest = max(math.dist(a, b), math.dist(b, c), math.dist(c, a))
    u = (b[0] - a[0], b[1] - a[1], b[2] - a[2])
    v = (c[0] - a[0], c[1] - a[1], c[2] - a[2])
    nx = u[1] * v[2] - u[2] * v[1]
    ny = u[2] * v[0] - u[0] * v[2]
    nz = u[0] * v[1] - u[1] * v[0]
    area2 = math.sqrt(nx * nx + ny * ny + nz * nz)
    return longest * longest / area2


def sides(tri):
    a, b, c = tri
    return math.dist(a, b), math.dist(b, c), math.dist(c, a)


@pytest.mark.parametrize("name", sorted(TRIS))
def test_conversion_is_the_closed_form(name):
    """library = legacy * p / (2*sqrt(3)*l_max), to round-off, on every triangle."""
    tri = TRIS[name]
    a, b, c = sides(tri)
    p = a + b + c
    lmax = max(a, b, c)
    predicted = legacy_aspect(tri) * p / (2.0 * math.sqrt(3.0) * lmax)
    assert triangle_quality(*tri).aspect_ratio == pytest.approx(predicted, rel=1e-14)


@pytest.mark.parametrize("name", sorted(TRIS))
def test_dialects_are_not_the_same_quantity(name):
    """The factor is bounded away from 1 on every triangle: the measures differ."""
    tri = TRIS[name]
    a, b, c = sides(tri)
    factor = (a + b + c) / (2.0 * math.sqrt(3.0) * max(a, b, c))
    assert 1.0 / math.sqrt(3.0) - 1e-12 <= factor <= math.sqrt(3.0) / 2.0 + 1e-12
    assert factor < 1.0 - 1e-3


def test_a_threshold_of_25_does_not_transfer():
    """Concrete counterexample: one triangle, two gates, opposite verdicts.

    ``(0,0)-(2,0)-(1,0.0741)``: legacy aspect 26.99 > 25 (the shipped
    ``mesh_conditioning`` gate rejects it) while the library aspect is 15.60 < 25
    (the new default gate accepts it).  Both gates use min angle 2 deg and this
    element's min angle is 4.24 deg, so the disagreement is purely the aspect
    dialect -- not a threshold that was tuned to manufacture a difference.
    """
    tri = ((0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (1.0, 0.0741, 0.0))
    legacy = legacy_aspect(tri)
    q = triangle_quality(*tri)
    assert legacy > LEGACY_LIMIT                  # old gate: fail
    assert q.aspect_ratio < LIBRARY_LIMIT          # new gate: pass
    assert q.min_angle_deg > 2.0                   # same angle threshold both sides
    assert q.aspect_ratio / legacy == pytest.approx(0.578, abs=2e-3)  # ~1/sqrt(3)
    # scaling cannot reconcile them: the gap is structural, not a geometric accident
    scaled = tuple((3.7 * x, 3.7 * y, 3.7 * z) for x, y, z in tri)
    assert legacy_aspect(scaled) == pytest.approx(legacy, rel=1e-12)
    assert triangle_quality(*scaled).aspect_ratio == pytest.approx(q.aspect_ratio, rel=1e-12)


def test_documented_gap_is_measured_in_the_artifact():
    """The audit artifact must exist and record the measured divergence, so the
    claim is not asserted in prose only."""
    import json
    from pathlib import Path

    path = Path(__file__).resolve().parent.parent / "artifacts" / "mesh_quality_audit" / "dialect.json"
    if not path.exists():
        pytest.skip("audit artifact not generated in this checkout")
    data = json.loads(path.read_text())
    assert data["verdict"]["same_quantity"] is False
    lo, hi = data["verdict"]["closed_form_ratio_range"]
    assert 0.577 < lo < hi <= 0.8661
    coarse = data["mesh_measurements"]["suboff_coarse"]
    assert coarse["elements_where_verdict_differs_at_limit_25"] == coarse[
        "elements_failing_legacy_but_passing_new"
    ]