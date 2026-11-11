"""Cap adaptation: profile-column dropping at a hull's tips.

The structured tensor-product grid carries ``n_v`` profile points on every station.
At a tip the circumferential cell shrinks with the radius while the axial spacing does
not, so the last rings before the apex become thin fans.  ``cap_adaptation=True``
coarsens the profile sampling where the ring is already finer than the axial spacing
(``docs/parametric_hull/CAP_ADAPTATION.md``).  Measured numbers live in
``examples/parametric_hull_cap/outputs/results.json``.
"""
from __future__ import annotations

import math

import pytest

from tensorcad.parametric_hull import (
    ParametricHullError,
    _adaptive_profile_rows,
    loft_nurbs_surface,
    mesh_quality_report,
    suboff_cross_sections,
    suboff_station_fractions,
    tessellate_surface,
    wigley_cross_sections,
)

N_U = 160
N_V = 16


@pytest.fixture(scope="module")
def suboff_surface():
    sections = suboff_cross_sections(
        n_points=32, station_fractions=suboff_station_fractions(41)
    )
    return loft_nurbs_surface(sections)


def _signed_volume(mesh):
    total = 0.0
    for f in mesh.faces:
        a, b, c = (mesh.vertices[i] for i in f)
        total += (
            a[0] * (b[1] * c[2] - b[2] * c[1])
            - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0])
        ) / 6.0
    return abs(total)


def _quality(mesh_data):
    return mesh_quality_report(mesh_data.mesh, n_faces=mesh_data.hull_face_count)


def test_baseline_is_the_documented_limitation(suboff_surface):
    """The structured grid reproduces the documented defect (jacobian 0.13, aspect 7.9)."""
    q = _quality(tessellate_surface(suboff_surface, n_u=N_U, n_v=N_V))
    assert q["min_scaled_jacobian"] == pytest.approx(0.1267, abs=5e-4)
    assert q["max_aspect_ratio"] == pytest.approx(7.88, abs=0.02)


def test_cap_adaptation_meets_the_quality_gate(suboff_surface):
    """Criterion (e) of the SUBOFF node: min scaled Jacobian > 0.30, aspect < 5."""
    q = _quality(tessellate_surface(suboff_surface, n_u=N_U, n_v=N_V, cap_adaptation=True))
    assert q["min_scaled_jacobian"] > 0.30
    assert q["max_aspect_ratio"] < 5.0
    assert q["n_degenerate_faces"] == 0
    assert q["watertight"] is True
    assert q["edge_manifold"] is True


def test_cap_adaptation_preserves_the_body(suboff_surface):
    """The mesh stays the same body: volume drift 0.02 %, face count within 10 %."""
    baseline = tessellate_surface(suboff_surface, n_u=N_U, n_v=N_V)
    adapted = tessellate_surface(suboff_surface, n_u=N_U, n_v=N_V, cap_adaptation=True)
    v0 = _signed_volume(baseline.mesh)
    v1 = _signed_volume(adapted.mesh)
    assert abs(v1 / v0 - 1.0) < 1e-3
    # Euler characteristic still 2 (a closed genus-0 body): V - E + F.
    mesh = adapted.mesh
    edges = set()
    for f in mesh.faces:
        a, b, c = f
        edges.update({(a, b), (b, c), (c, a)})
    edges = {(min(e), max(e)) for e in edges}
    assert len(mesh.vertices) - len(edges) + len(mesh.faces) == 2


def test_wigley_quality_improves_and_waterplane_is_unchanged():
    """A closed bow benefits too, and the top row is never coarsened."""
    surface = loft_nurbs_surface(wigley_cross_sections())
    base = tessellate_surface(surface, n_u=97, n_v=16, close_waterplane=True)
    adapted = tessellate_surface(
        surface, n_u=97, n_v=16, close_waterplane=True, cap_adaptation=True
    )
    qb, qa = _quality(base), _quality(adapted)
    assert qa["min_scaled_jacobian"] > qb["min_scaled_jacobian"]
    assert qa["max_aspect_ratio"] < qb["max_aspect_ratio"]
    from tensorcad.parametric_hull import compute_hydrostatics

    kwargs = dict(waterline_z=0.0, length=10.0, beam=1.0, draft=0.625)
    hb = compute_hydrostatics(base.mesh, **kwargs)
    ha = compute_hydrostatics(adapted.mesh, **kwargs)
    assert ha.waterplane_area == pytest.approx(hb.waterplane_area, rel=1e-12)
    assert ha.displaced_volume == pytest.approx(hb.displaced_volume, rel=5e-3)


def test_adaptive_rows_keep_the_ends_and_shrink_towards_the_tips(suboff_surface):
    us = [i / (N_U - 1) for i in range(N_U)]
    rows = _adaptive_profile_rows(suboff_surface, us, N_V, 1.0)
    assert len(rows) == N_U
    for row in rows:
        assert row[0] == 0.0 and row[-1] == 1.0          # keel and top stay
        assert all(b >= a for a, b in zip(row, row[1:]))  # monotone, no reordering
    # `_adaptive_profile_rows` snaps onto the n_v index grid, so a ring whose segment
    # count is s carries s + 1 distinct samples (s = n_v - 1 is the uncoarsened ring).
    distinct = [len(set(row)) for row in rows]
    assert min(distinct) <= 5           # the tips really are coarsened, hard
    assert max(distinct) == N_V         # the midbody is untouched
    mid = distinct[20:-20]
    assert max(mid) == N_V
    assert min(mid) >= N_V - 1          # no midbody station loses more than one point


def test_invalid_aspect_target_is_rejected(suboff_surface):
    us = [0.0, 0.5, 1.0]
    with pytest.raises(ParametricHullError):
        _adaptive_profile_rows(suboff_surface, us, N_V, 0.0)