"""Contracts of ``tensorcad.geomquality.mesh_defects`` on meshes whose answer is known.

No corpus is used here on purpose: the corpus claims live in
``examples/mesh_route_naming`` (they need the BRep injections). These tests build the
smallest mesh that carries each defect **sewn into one component**, so the codes cannot be
earned through the ``ISOLATED_PATCH`` shortcut -- a sliver or a micro edge in a real
re-tessellation is part of the shell, not a detached island.
"""
from __future__ import annotations

import pytest

from tensorcad.geometry import Mesh
from tensorcad.geomquality import (
    name_mesh_defects,
    repair_facet_defects,
)


def tetrahedron(scale: float = 1.0) -> Mesh:
    """A closed, outward-oriented unit tetrahedron -- four faces, six edges."""
    vertices = [
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0),
    ]
    faces = [(0, 2, 1), (0, 1, 3), (0, 3, 2), (1, 2, 3)]
    scaled = [(scale * x, scale * y, scale * z) for x, y, z in vertices]
    return Mesh(vertices=tuple(scaled), faces=tuple(faces))


def sewn_sliver() -> Mesh:
    """One ordinary triangle plus a 1e-9 m high one sharing its edge: one component."""
    vertices = [
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (0.5, 1e-9, 0.0),  # the needle vertex
        (0.5, 0.0, 1.0),
    ]
    faces = [(0, 1, 2), (0, 1, 3)]
    return Mesh(vertices=tuple(vertices), faces=tuple(faces))


def sewn_micro_edge() -> Mesh:
    """A 1e-9 x 1 m band: every edge of the short side is below the weld tolerance."""
    vertices = [
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (1.0, 1e-9, 0.0),
        (0.0, 1e-9, 0.0),
    ]
    faces = [(0, 1, 2), (0, 2, 3)]
    return Mesh(vertices=tuple(vertices), faces=tuple(faces))


def with_detached_triangle(mesh: Mesh, height: float = 1.0) -> Mesh:
    """``mesh`` plus an island that shares no edge with it.

    ``height`` sets the island's area; the default island is an *ordinary* triangle, so it
    survives every facet-scale criterion and is removed by nothing the repair does.
    """
    offset = len(mesh.vertices)
    vertices = list(mesh.vertices) + [
        (5.0, 0.0, 0.0),
        (6.0, 0.0, 0.0),
        (5.0, height, 0.0),
    ]
    faces = list(mesh.faces) + [(offset, offset + 1, offset + 2)]
    return Mesh(vertices=tuple(vertices), faces=tuple(faces))


# --------------------------------------------------------------------------- #
# the clean reference
# --------------------------------------------------------------------------- #
def test_a_closed_tetrahedron_has_no_codes() -> None:
    report = name_mesh_defects(tetrahedron(), declared_length=1.0)
    assert report.codes == ()
    assert report.clean
    assert report.metrics["faces"] == 4
    assert report.metrics["open_edges"] == 0
    assert report.metrics["components"] == [4]


def test_open_and_nonmanifold_edges_are_context_not_codes() -> None:
    """Topology belongs to ``audit_watertightness``; this module must not re-dialect it."""
    report = name_mesh_defects(sewn_sliver(), declared_length=1.0)
    assert report.metrics["open_edges"] == 4
    assert "FREE_EDGE" not in report.codes
    assert "OPEN_TOP" not in report.codes


# --------------------------------------------------------------------------- #
# the two facet-scale codes, earned without the island shortcut
# --------------------------------------------------------------------------- #
def test_a_sewn_sliver_is_named_by_its_area_alone() -> None:
    mesh = sewn_sliver()
    report = name_mesh_defects(mesh, declared_length=1.0)
    assert report.codes == ("SLIVER",)
    assert report.metrics["components"] == [2]  # one component: no island to lean on
    defect = report.defect("SLIVER")
    assert defect is not None
    assert defect.count == 1
    assert defect.measure == pytest.approx(0.5)
    assert defect.limit == pytest.approx(0.02)
    assert "m^2" in defect.detail


def test_a_sewn_micro_edge_is_named_by_its_length_alone() -> None:
    report = name_mesh_defects(sewn_micro_edge(), declared_length=1.0, sliver_share_limit=1.0)
    assert "MICRO_EDGE" in report.codes
    assert report.metrics["components"] == [2]
    defect = report.defect("MICRO_EDGE")
    assert defect is not None
    assert defect.count == 2
    assert defect.measure == pytest.approx(0.4)
    assert defect.measure > defect.limit


def test_the_needle_is_named_under_both_limits_that_it_earns() -> None:
    """A 1e-9 x 1 m band is a sliver *and* carries sub-tolerance edges; both are reported."""
    report = name_mesh_defects(sewn_micro_edge(), declared_length=1.0)
    assert report.codes == ("MICRO_EDGE", "SLIVER")


def test_an_island_is_named_isolated_patch() -> None:
    report = name_mesh_defects(with_detached_triangle(tetrahedron()), declared_length=1.0)
    assert report.codes == ("ISOLATED_PATCH",)
    defect = report.defect("ISOLATED_PATCH")
    assert defect is not None
    assert defect.count == 1  # one minor component
    assert report.metrics["components"] == [4, 1]


# --------------------------------------------------------------------------- #
# limits: where the criteria stop firing
# --------------------------------------------------------------------------- #
def test_a_sliver_share_at_the_limit_does_not_fire() -> None:
    """Two faces, one sliver: 50 % fires; raise the limit to 50 % and the code is gone."""
    assert name_mesh_defects(sewn_sliver(), declared_length=1.0).named("SLIVER")
    quiet = name_mesh_defects(sewn_sliver(), declared_length=1.0, sliver_share_limit=0.5)
    assert not quiet.named("SLIVER")


def test_without_a_declared_length_the_sliver_criterion_is_not_evaluated() -> None:
    """No scale, no sliver limit: the report says so instead of inventing one."""
    report = name_mesh_defects(sewn_sliver())
    assert report.sliver_area_limit is None
    assert "SLIVER" not in report.codes
    assert report.metrics["sliver_share"] == 0.0


def test_the_share_criteria_are_scale_free_but_the_weld_tolerance_is_physical() -> None:
    """Scaling the model 1000x keeps the codes only if the tolerance is scaled with it."""
    base = name_mesh_defects(sewn_micro_edge(), declared_length=1.0)
    scaled = name_mesh_defects(
        sewn_micro_edge(), declared_length=1.0
    )  # same object, for the identity check
    assert base.codes == scaled.codes
    big = Mesh(
        vertices=tuple(tuple(1000.0 * c for c in point) for point in sewn_micro_edge().vertices),
        faces=sewn_micro_edge().faces,
    )
    assert "MICRO_EDGE" not in name_mesh_defects(big, declared_length=1000.0).codes
    assert "MICRO_EDGE" in name_mesh_defects(
        big, declared_length=1000.0, weld_tolerance=1e-3
    ).codes
    assert "SLIVER" in name_mesh_defects(big, declared_length=1000.0).codes


def test_bad_inputs_are_rejected() -> None:
    with pytest.raises(ValueError):
        name_mesh_defects(Mesh(vertices=((0.0, 0.0, 0.0),), faces=()))
    with pytest.raises(ValueError):
        name_mesh_defects(tetrahedron(), declared_length=1.0, weld_tolerance=0.0)
    with pytest.raises(ValueError):
        name_mesh_defects(tetrahedron(), declared_length=-1.0)


# --------------------------------------------------------------------------- #
# repair: it either improves the mesh or says why it did not
# --------------------------------------------------------------------------- #
def test_repair_of_a_clean_mesh_is_a_no_op() -> None:
    result = repair_facet_defects(tetrahedron(), declared_length=1.0)
    assert result.accepted
    assert result.codes_before == ()
    assert result.codes_after == ()
    assert result.merged_vertices == 0
    assert result.removed_faces == 0
    assert result.volume_relative_change == pytest.approx(0.0)
    assert result.area_relative_change == pytest.approx(0.0)


def test_repair_clears_a_degenerate_island_and_leaves_the_shell_alone() -> None:
    """Safe case: the island is itself below the sliver area, so dropping it costs nothing."""
    mesh = with_detached_triangle(tetrahedron(), height=1e-10)
    result = repair_facet_defects(mesh, declared_length=1.0)
    assert set(result.codes_before) == {"ISOLATED_PATCH", "MICRO_EDGE", "SLIVER"}
    assert result.codes_after == ()
    assert result.accepted
    assert result.removed_faces == 1
    # the island carries no volume, so the displacement must not move
    assert result.volume_relative_change < 1e-9
    assert result.area_relative_change < 1e-9


def test_repair_of_the_sewn_needle_either_clears_it_or_refuses_honestly() -> None:
    """Collapsing a 1e-9 m band destroys both faces: the refusal is the honest outcome."""
    result = repair_facet_defects(sewn_micro_edge(), declared_length=1.0)
    if result.accepted:
        assert result.codes_after == ()
    else:
        assert result.mesh is not None
        assert "original mesh returned" in result.reason
        assert result.codes_after == result.codes_before


def test_repair_refuses_rather_than_returning_nothing() -> None:
    """A tolerance that eats the model is a refusal, and the original mesh comes back."""
    original = tetrahedron()
    result = repair_facet_defects(original, declared_length=1.0, weld_tolerance=10.0)
    assert not result.accepted
    assert "every face" in result.reason
    assert result.mesh is original


def test_repair_refuses_when_a_measured_quantity_would_move() -> None:
    """Whatever the guard trips on, it is named: a silent change of shape is not allowed."""
    result = repair_facet_defects(
        sewn_sliver(),
        declared_length=1.0,
        sliver_share_limit=0.0,
        volume_tolerance=0.0,
        area_tolerance=0.0,
    )
    assert not result.accepted
    assert ("volume moved" in result.reason) or ("area moved" in result.reason)
    assert result.reason.endswith("original mesh returned")


# --------------------------------------------------------------------------- #
# serialisation
# --------------------------------------------------------------------------- #
def test_the_report_serialises_with_its_evidence() -> None:
    report = name_mesh_defects(sewn_sliver(), declared_length=1.0)
    payload = report.to_dict()
    assert payload["codes"] == ["SLIVER"]
    assert payload["defects"][0]["code"] == "SLIVER"
    assert payload["defects"][0]["repair"]
    assert payload["metrics"]["faces"] == 2
    assert payload["declared_length"] == 1.0
    for key in ("weld_tolerance", "sliver_area_limit", "sliver_share_limit"):
        assert key in payload


def test_an_ordinary_island_is_reported_and_not_silently_deleted() -> None:
    """Welding cannot sew an island, so the repair refuses instead of deleting geometry."""
    result = repair_facet_defects(with_detached_triangle(tetrahedron()), declared_length=1.0)
    assert result.codes_before == ("ISOLATED_PATCH",)
    assert result.codes_after == ("ISOLATED_PATCH",)
    assert not result.accepted
    assert result.mesh is not None
    assert "no facet-scale repair available for ISOLATED_PATCH" in result.reason


def test_the_repair_result_serialises_with_its_outcome() -> None:
    result = repair_facet_defects(
        with_detached_triangle(tetrahedron(), height=1e-10), declared_length=1.0
    )
    payload = result.to_dict()
    assert payload["codes_after"] == []
    assert payload["accepted"] is True
    assert "welded" in payload["reason"]


def test_the_module_needs_no_third_party_stack() -> None:
    """The package's zero-dependency claim has to survive this module."""
    module_source = __import__("tensorcad.geomquality.mesh_defects", fromlist=["x"])
    imports = {
        line.split()[1].split(".")[0]
        for line in open(module_source.__file__, encoding="utf-8")
        if line.startswith("import ") or line.startswith("from ")
    }
    assert imports <= {"__future__", "math", "collections", "dataclasses", "typing", "..", ""}
    assert "numpy" not in imports and "torch" not in imports
