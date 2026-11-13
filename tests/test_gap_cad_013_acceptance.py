"""Machine gate for the ``GAP-CAD-013`` acceptance criterion.

The anchor in ``docs/oss_inventory_cad/GAP-ANCHORS.md`` states the criterion as a
sentence that a machine can adjudicate:

    中面抽取：对厚度 t 的板件抽取中面，面积误差 <1%，面数 ≤ 原面数 1/10；
    特征抑制后网格规模下降 ≥50%。

This file turns each of the three clauses into one measurement and asserts the
threshold *the anchor states*, not the value observed on the day the anchor was
written.  ``tests/test_midsurface.py`` already pins the physics (analytic mid-plane
areas, curvature surplus, hull plating); this file pins the **acceptance contract**
of the gap ledger, so flipping the ledger status away from ``absent`` is backed by a
test that fails if the contract regresses.

Clause 2 is measured at the stage that actually delivers it.  Mid-surface extraction
keeps one facet per wall -- a plate wall is two paired source faces -- so it halves the
face count and no more (asserted below); the ``≤ 1/10`` face budget is therefore met by
the planar-defeaturing stage of the same feature --
which is why the ledger records the feature as ``partial`` with that boundary spelled
out, rather than claiming extraction alone deduplicates the geometry.
"""
from __future__ import annotations

import pytest

from tensorcad import mesh as meshlib
from tensorcad.midsurface import (
    defeature_planar_regions,
    extract_midsurface,
    plate_solid,
    cylinder_skin,
)

LENGTH, BREADTH, THICKNESS = 8.0, 3.0, 0.02
PLATE = plate_solid(LENGTH, BREADTH, THICKNESS, nx=20, ny=10)

#: Clause thresholds, quoted from the anchor (not from hardware behaviour).
AREA_TOLERANCE = 0.01  # <1 % area error
FACE_BUDGET = 1.0 / 10.0  # 面数 ≤ 原面数 1/10
SCALE_REDUCTION = 0.50  # 网格规模下降 ≥50 %

MID_SURFACE = extract_midsurface(PLATE)
IDEALISED = defeature_planar_regions(MID_SURFACE.mesh)


# --------------------------------------------------------------------------- #
# clause 1 -- area error < 1 % on a plate of thickness t
# --------------------------------------------------------------------------- #
def test_clause1_mid_surface_area_is_within_one_percent():
    """The mid-plane of a flat wall is the mid-plane, so the reference is exact."""
    reference = LENGTH * BREADTH
    error = (MID_SURFACE.area - reference) / reference
    assert abs(error) < AREA_TOLERANCE, f"area error {error:+.3e} exceeds {AREA_TOLERANCE:.0%}"
    # and the extraction measured the wall it was given, rather than assuming t
    assert MID_SURFACE.mean_thickness == pytest.approx(THICKNESS, rel=1e-9)


# --------------------------------------------------------------------------- #
# clause 2 -- face count <= 1/10 of the source
# --------------------------------------------------------------------------- #
def test_clause2_extraction_alone_keeps_one_facet_per_paired_face():
    """Boundary of the claim: extraction does not deduplicate geometry.

    Every paired source face contributes exactly one mid-surface facet, so clause 2
    cannot be satisfied by the extraction stage on its own -- the ledger must not
    pretend otherwise.
    """
    # one facet per *wall* -- a plate wall has two paired source faces
    assert 2 * len(MID_SURFACE.mesh.faces) == MID_SURFACE.paired_face_count
    assert MID_SURFACE.paired_face_count + len(MID_SURFACE.unpaired) == MID_SURFACE.source_face_count
    assert len(MID_SURFACE.mesh.faces) > FACE_BUDGET * MID_SURFACE.source_face_count


def test_clause2_idealised_face_count_is_within_the_face_budget():
    """After planar defeaturing the plate is two triangles: 2 / 400 = 0.5 %."""
    source_faces = MID_SURFACE.source_face_count
    idealised_faces = IDEALISED.face_count_after
    ratio = idealised_faces / source_faces
    assert ratio <= FACE_BUDGET, (
        f"{idealised_faces} faces vs {source_faces} source faces = {ratio:.3%}, "
        f"budget {FACE_BUDGET:.0%}"
    )
    assert idealised_faces == 2  # the two sides of one flat wall
    assert len(IDEALISED.regions) == 1


# --------------------------------------------------------------------------- #
# clause 3 -- mesh scale down by >= 50 % after feature suppression
# --------------------------------------------------------------------------- #
def test_clause3_face_count_reduction_exceeds_fifty_percent():
    assert IDEALISED.face_reduction >= SCALE_REDUCTION, (
        f"face reduction {IDEALISED.face_reduction:.3%} < {SCALE_REDUCTION:.0%}"
    )


def test_clause3_vertex_count_reduction_exceeds_fifty_percent():
    before = len(MID_SURFACE.mesh.vertices)
    after = len(IDEALISED.mesh.vertices)
    reduction = 1.0 - after / before
    assert reduction >= SCALE_REDUCTION, (
        f"vertex reduction {reduction:.3%} < {SCALE_REDUCTION:.0%}"
    )


# --------------------------------------------------------------------------- #
# the suppression must not buy the reduction with geometry
# --------------------------------------------------------------------------- #
def test_defeature_preserves_area_and_is_watertight_free():
    """Face budget is spent on coplanar merging only: area is preserved to round-off."""
    assert IDEALISED.area_relative_error < 1e-12
    assert IDEALISED.mesh.faces
    assert not meshlib.is_watertight(IDEALISED.mesh)  # a mid-surface is open by nature


def test_curved_skin_is_not_suppressed_and_never_drifted():
    """A curved skin has no coplanar neighbours: no reduction is claimed, none invented."""
    skin = cylinder_skin(radius=1.0, thickness=0.02, length=4.0, segments=64)
    mid = extract_midsurface(skin)
    result = defeature_planar_regions(mid.mesh)
    assert result.face_reduction == 0.0
    # the only merger allowed is a strip of two coplanar facets along the axis
    assert max(len(region.faces) for region in result.regions) == 2
    assert all(len(region.boundary_loops) == 1 for region in result.regions)
    assert result.area_relative_error < 1e-12

# --------------------------------------------------------------------------- #
# the ledger and the code may not drift apart
# --------------------------------------------------------------------------- #
def test_ledger_records_the_landed_capability():
    """``feature_points.json`` + the rendered anchor must agree, and cite live files.

    The anchor is rendered by ``tools/oss_inventory_cad/build.py``; this gate keeps the
    claim ("partial, here is the evidence") tied to files that actually exist, so the
    ledger cannot drift back to ``absent`` while the code is still there -- or claim
    ``real`` without the acceptance clause above being met.
    """
    import json
    import pathlib

    root = pathlib.Path(__file__).resolve().parents[1]
    points = json.loads((root / "tools/oss_inventory_cad/feature_points.json").read_text())
    f13 = next(f for f in points["feature_points"] if f["id"] == "F13")
    record = f13["tensorcad"]

    assert record["status"] == "partial", "F13 landed; it is partial, not absent and not real"
    assert "tensorcad/midsurface.py" in record["evidence"]
    for evidence in record["evidence"]:
        assert (root / evidence).exists(), f"cited evidence is missing: {evidence}"

    anchor = (root / "docs/oss_inventory_cad/GAP-ANCHORS.md").read_text()
    block = anchor[anchor.index("### GAP-CAD-013"):anchor.index("### GAP-CAD-014")]
    assert "`partial`" in block
    assert record["note"] in block, "rendered anchor drifted from feature_points.json"
    for evidence in record["evidence"]:
        assert f"`{evidence}`" in block
