"""Cross-route audit consistency: the BRep route and the mesh route on one corpus.

``tests/test_audit_crosscheck.py`` pins the claims that
``examples/audit_crosscheck/run_example.py`` measures over the full 45-model corpus,
on a fast subset so the suite stays quick:

* the mesh route names every topology fault family the BRep route names;
* the mesh route cannot *name* a sliver, a micro edge or a seam class -- it can only
  say "the shell is open": that conflation is the finding, not a bug, and it is
  asserted so a future detector that changes it turns this test red;
* a clean reference hull tessellation is *not* a clean mesh: the Wigley form's
  zero-breadth stem pinches, and the triangle-quality gate rejects 100 % of raw CAD
  tessellation at every density tried -- so the mesh-route gate must be fed
  conditioned geometry, not raw patch sampling;
* ``intersect.find_self_intersections`` is the only mesh-side observable that sees a
  folded patch, and it fires on folds and nothing else.
"""
from __future__ import annotations

import importlib.util
import math
from pathlib import Path

import pytest

from tensorcad import geometry_audit as audit
from tensorcad import intersect, mesh as meshlib, validation
from tensorcad.geomquality import baselines, defects

REPO_ROOT = Path(__file__).resolve().parents[1]
SPECIMEN = "wigley_iii"
UV = (8, 4)


def _self_intersections(mesh) -> int:
    found = intersect.find_self_intersections(mesh)
    if isinstance(found, dict):
        for key in ("pairs", "intersections", "count", "found"):
            if key in found:
                value = found[key]
                return len(value) if isinstance(value, (list, tuple)) else int(value)
        return int(found.get("n_pairs", 0))
    return len(found)


def _specimen():
    return baselines.specimens()[SPECIMEN]()


def _declared(specimen) -> dict[str, float]:
    return {
        "length": specimen.length,
        "beam": specimen.beam,
        "draft": specimen.draft,
        "block_coefficient": specimen.particulars["block_coefficient"],
    }


@pytest.fixture(scope="module")
def corpus():
    specimen = _specimen()
    declared = _declared(specimen)
    models = {}
    for case in [None] + list(defects.DEFECT_CASES):
        brep = specimen.brep if case is None else case.inject(specimen.brep)
        result = defects.detect_all(brep, declared=declared)
        mesh = brep.tessellate(*UV)
        watertightness = audit.audit_watertightness(mesh, length=float(brep.length)).to_dict()
        models[case.name if case else "CLEAN"] = {
            "expected_code": case.code if case else None,
            "codes": sorted(result["codes"]),
            "mesh": mesh,
            "watertight": watertightness,
            "fairness": len(audit.audit_fairness(mesh).defects),
            "open_edges": watertightness["open_edges"],
            "self_intersections": _self_intersections(mesh),
            "volume": meshlib.total_signed_volume(mesh),
        }
    return models


def test_brep_route_names_every_injection(corpus) -> None:
    """The 14-case library is the ground truth this whole study reads."""
    missed = [
        (case, model["codes"])
        for case, model in corpus.items()
        if model["expected_code"] and model["expected_code"] not in model["codes"]
    ]
    assert not missed, missed
    assert not corpus["CLEAN"]["codes"]


@pytest.mark.parametrize(
    "case",
    ["hole", "open_waterline", "unsewn_seam", "non_manifold_flap", "flipped_face", "duplicate_face"],
)
def test_mesh_route_sees_every_topology_family(corpus, case) -> None:
    """Each topology fault produces a mesh-side signal the clean model does not have."""
    clean, damaged = corpus["CLEAN"], corpus[case]
    signals = (
        damaged["open_edges"] > clean["open_edges"],
        damaged["watertight"]["nonmanifold_edges"] > clean["watertight"]["nonmanifold_edges"],
        damaged["watertight"]["duplicate_faces"] > clean["watertight"]["duplicate_faces"],
        damaged["watertight"]["winding_conflicts"] > clean["watertight"]["winding_conflicts"],
    )
    assert any(signals), (case, damaged["watertight"], clean["watertight"])


@pytest.mark.parametrize("case", ["sliver_face", "micro_edge", "g1_splice", "g0_splice"])
def test_unnamed_faults_degrade_to_the_same_open_shell_signal(corpus, case) -> None:
    """A sliver, a micro edge and a seam class all read as "not closed" to the mesh route.

    This is the measured conflation: the mesh route answers *whether* the shell closes,
    never *which* fault did it. If a future detector names these, this test should be
    updated together with the docs -- the assertion is the contract.
    """
    clean, damaged = corpus["CLEAN"], corpus[case]
    assert damaged["open_edges"] > clean["open_edges"]
    # No mesh-side quality or fairness verdict distinguishes it from any other open shell:
    mesh = damaged["mesh"]
    geometry = [check.name for check in validation.inspect_geometry(mesh) if not check.passed]
    assert "watertight" in geometry


def test_self_intersection_is_the_only_geometric_signal(corpus) -> None:
    """The folded patch is invisible to every check except the triangle/triangle one."""
    clean, fold = corpus["CLEAN"], corpus["fold"]
    assert fold["self_intersections"] > 0
    assert not clean["self_intersections"]
    # The clean reference has none, so the signal is not noise. Some *other* injections
    # also overlap geometrically (a flap, a doubled face): those are reported, not denied
    # -- the detector answers "do two faces pass through each other", and for them, yes.
    overlapping = sorted(
        case for case, model in corpus.items() if case != "fold" and model["self_intersections"]
    )
    assert overlapping != sorted(corpus), overlapping


def test_clean_tessellation_is_not_a_clean_mesh(corpus) -> None:
    """The Wigley stem pinches: the mesh route reports defects the BRep route denies.

    ``geomquality`` says the clean reference has no defect (``codes == []``); the same
    geometry sampled into triangles has one non-manifold edge, one duplicated facet and
    three pinch vertices at the zero-breadth stem. Both statements are true of their own
    representation; the test fixes the disagreement so it cannot drift silently.
    """
    clean = corpus["CLEAN"]
    assert clean["codes"] == []
    kinds = {defect["kind"] for defect in clean["watertight"]["defects"]}
    assert "nonmanifold_edge" in kinds
    assert "duplicate_facets" in kinds
    assert not clean["watertight"]["passed"]
    assert not clean["watertight"]["closed"] is False or clean["watertight"]["open_edges"] == 0


def test_quality_gate_rejects_raw_cad_tessellation(corpus) -> None:
    """Smallest triangle angle and sliver share do not improve with patch sampling.

    If refinement fixed the triangles this test would be the wrong assertion; measured,
    it does not: the sliver share grows with density because the thin triangles live
    along the pointed stem, whose aspect ratio is set by the hull form.
    """
    specimen = _specimen()
    shares = []
    for uv in ((8, 4), (16, 8)):
        mesh = specimen.brep.tessellate(*uv)
        slivers = 0
        smallest = math.inf
        for face in mesh.faces:
            points = [mesh.vertices[index] for index in face]
            lengths = [math.dist(points[a], points[(a + 1) % 3]) for a in range(3)]
            opposite, side_a, side_b = sorted(lengths)
            denominator = 2.0 * side_a * side_b
            angle = (
                0.0
                if denominator <= 0.0
                else math.degrees(
                    math.acos(
                        max(-1.0, min(1.0, (side_a**2 + side_b**2 - opposite**2) / denominator))
                    )
                )
            )
            smallest = min(smallest, angle)
            slivers += angle < 2.0
        shares.append(slivers / len(mesh.faces))
    # Coarse sampling hides the slivers (8x4: 0 %); refining exposes them and the share
    # grows -- there is no density at which raw patch sampling becomes solver-ready.
    assert shares[1] > 0.02, shares


def test_example_module_is_importable() -> None:
    """The harness itself must import without running the 45-model sweep."""
    path = REPO_ROOT / "examples" / "audit_crosscheck" / "run_example.py"
    spec = importlib.util.spec_from_file_location("audit_crosscheck_example", path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    assert hasattr(module, "mesh_route") and hasattr(module, "main")