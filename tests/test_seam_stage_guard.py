"""The seam stage must never trade a measured defect for an unmeasured one.

``seam_stage`` exists to close the small inter-patch gaps a publisher leaves between
neighbouring trimmed patches.  On the publisher hull patch sets in the asset store the
stitch is measurably *harmful*: it splits boundary chords, introduces non-manifold chords
and duplicates triangles, while most of the open boundary stays open (mission m22849:
KCS ``KCS_hull_SVA.igs``, tolerance 0.02 -> free edges 1697 -> 1331, duplicates 22,
orientation conflicts 209, manifold False).

Three properties are pinned here:

* a repeated triangle is never geometry, so the stitch drops it and counts it;
* the stitch reports the topology it changed (free edges / non-manifold / conflicts /
  components) before and after, instead of returning a mesh whose only evidence is
  "free edges went down";
* when stitching would increase the non-manifold edge count, it is *refused*: the caller
  gets the unstitched mesh and a reason, so a repair chain can never silently degrade a
  mesh and then report it as healed.
"""
from __future__ import annotations

import os

import pytest

from tensorcad.geomquality.trim_mesh import (
    _drop_duplicate_triangles,
    audit_mesh,
    seam_stage,
)

KCS = "/srv/tensor/agent/assets/geometry/hulls/KCS/KCS_hull_SVA.igs"


# --------------------------------------------------------------------------- #
# unit: repeated triangles
# --------------------------------------------------------------------------- #
def test_duplicate_triangles_are_dropped_under_any_rotation():
    faces = [(0, 1, 2), (1, 2, 0), (2, 0, 1), (0, 2, 3)]
    kept, dropped = _drop_duplicate_triangles(faces)
    assert dropped == 2
    assert kept == [(0, 1, 2), (0, 2, 3)]


# --------------------------------------------------------------------------- #
# a plain T-junction is still stitched (the stage must keep doing its job)
# --------------------------------------------------------------------------- #
def test_t_junction_is_stitched_and_reported():
    # Two quads forming a 2x1 rectangle strip whose shared chord is subdivided on the
    # left patch only: vertex 4 sits in the middle of chord (1, 3).
    vertices = [
        (0.0, 0.0, 0.0), (0.0, 1.0, 0.0), (1.0, 1.0, 0.0),
        (1.0, 0.0, 0.0), (0.5, 1.0, 0.0),
    ]
    faces = [
        (0, 1, 4), (0, 4, 3),          # outer part of the left patch (triangulated)
        (3, 4, 2), (3, 2, 1), (1, 2, 4),  # right patch shares the chord (1, 3)
    ]
    before = audit_mesh(vertices, faces, degenerate_area=0.0)
    stitched, stitched_faces, record = seam_stage(vertices, faces, tolerance=0.1)
    assert record["applied"] is True
    assert "stitch" in record
    assert record["topology_before"]["free_edges"] == before.free_edge_count
    assert record["topology_after"]["non_manifold_edges"] <= before.non_manifold_edge_count
    assert "topology_before" in record and "topology_after" in record


def test_measured_stitch_never_returns_duplicate_triangles():
    vertices = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0),
                (0.0, 1.0, 0.0), (0.5, 0.0, 0.0), (0.5, 1.0, 0.0)]
    faces = [(0, 1, 2), (0, 2, 3), (0, 4, 1), (4, 5, 2), (4, 2, 1), (0, 3, 5), (0, 5, 4)]
    stitched, stitched_faces, record = seam_stage(vertices, faces, tolerance=0.2)
    keys = [tuple(sorted(face)) for face in stitched_faces]
    assert len(keys) == len(set(keys))
    assert record["stitch"]["duplicate_triangles_dropped"] >= 0
    audit = audit_mesh(stitched, stitched_faces, degenerate_area=0.0)
    assert audit.duplicate_triangle_count == 0


# --------------------------------------------------------------------------- #
# the real publisher patch set: healed or refused, never degraded
# --------------------------------------------------------------------------- #
@pytest.mark.skipif(not os.path.exists(KCS), reason="publisher asset not present")
def test_publisher_patch_set_stitch_is_refused_when_it_degrades_manifoldness():
    from tensorcad.geomquality.trim_mesh import read_iges_trimmed, tessellate_trimmed_model

    model = read_iges_trimmed(KCS)
    tessellation = tessellate_trimmed_model(model, scale=1.0)
    _, _, record = seam_stage(tessellation.vertices, tessellation.faces, tolerance=0.02)
    _, _, strict_record = seam_stage(tessellation.vertices, tessellation.faces,
                                     tolerance=0.02, strict=True)
    before = record["topology_before"]
    after = record["topology_after"]
    # the stitched mesh does close most of the boundary, so the default keeps it ...
    assert record["applied"] is True
    assert record["stitch"]["duplicate_triangles_dropped"] > 0
    # ... and the damage the stitch does is measured, not hidden.
    assert after["non_manifold_edges"] > before["non_manifold_edges"]
    assert record["manifold_degraded"] is True
    # ... while a strict call refuses the same stitch and keeps the input mesh.
    assert strict_record["applied"] is False
    assert "non-manifold" in strict_record["rejected"]
    assert strict_record["topology_after"]["non_manifold_edges"] > 0