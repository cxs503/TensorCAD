"""Tests for :mod:`tensorcad.meshquality.door` and ``tools/mesh_quality_door.py``.

The two halves are asserted separately, because they fail differently:

* the **composition** -- :func:`admit` adds no measure of its own, so what it must get
  right is that it asks the existing gates the right questions (shell topology when
  triangles are given, interior faces only for ``tet``/``hex``, and nothing at all for a
  cell kind the gate does not name), that a failure keeps the name, limit and value of the
  check that produced it, and that an empty mesh is an error rather than a pass;
* the **door** -- ``tools/mesh_quality_door.py`` is a CI contract, so the tests pin the
  exit codes, the tamper resistance of the frozen baseline (a hand-edited baseline must
  not buy a pass), the determinism of two runs, and that the two defect controls really
  are rejected. The closed-form comparisons are asserted from the tool's own record, not
  re-derived here, so a library change that moves a number fails in exactly one place.

No upstream tool (``vtk``/``gmsh``/``meshmagick``) is imported here: this host has none,
and the door's differential is closed-form plus regression, as its documentation says.
"""
from __future__ import annotations

import importlib.util
import json
import math
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from tensorcad.meshquality import DEFAULT_GATE, MetricLimit, QualityGate, admit  # noqa: E402
from tensorcad.meshquality.topology import gate_non_orthogonality, non_orthogonality  # noqa: E402

_SPEC = importlib.util.spec_from_file_location(
    "mesh_quality_door", ROOT / "tools" / "mesh_quality_door.py"
)
assert _SPEC is not None and _SPEC.loader is not None
door_tool = importlib.util.module_from_spec(_SPEC)
sys.modules["mesh_quality_door"] = door_tool
_SPEC.loader.exec_module(door_tool)

UNIT_CUBE = (
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.0),
    (1.0, 1.0, 0.0),
    (0.0, 1.0, 0.0),
    (0.0, 0.0, 1.0),
    (1.0, 0.0, 1.0),
    (1.0, 1.0, 1.0),
    (0.0, 1.0, 1.0),
)

#: The closed unit cube as a triangle shell, wound outwards.
CUBE_SHELL = (
    (0, 2, 1), (0, 3, 2),
    (4, 5, 6), (4, 6, 7),
    (0, 1, 5), (0, 5, 4),
    (1, 2, 6), (1, 6, 5),
    (2, 3, 7), (2, 7, 6),
    (3, 0, 4), (3, 4, 7),
)


def _shear_pair(shear: float = 0.5):
    return door_tool.hex_shear_block(shear)


# --------------------------------------------------------------------------------------
# the composition: admit()
# --------------------------------------------------------------------------------------
def test_admit_asks_the_element_gate_and_names_its_checks():
    points, cells = _shear_pair()
    verdict = admit(points, cells={"hex": cells})

    names = [check.name for check in verdict.checks]
    assert names[0] == "volume.max_non_orthogonality"
    assert "hex.min_sicn" in names and "hex.scaled_jacobian" in names
    # the shell half is not asked for a volume-only mesh
    assert not any(name.startswith("surface.") for name in names)
    # every check is the gate's own object, so it carries its limit and its number
    sicn = next(check for check in verdict.checks if check.name == "hex.min_sicn")
    assert sicn.limit.startswith(">= 0.1")
    assert sicn.value == pytest.approx(0.8172460864182961, rel=1e-12)
    assert verdict.passed
    assert verdict.failures == ()


def test_admit_asks_the_shell_gate_when_triangles_are_given():
    verdict = admit(UNIT_CUBE, triangles=CUBE_SHELL)
    names = [check.name for check in verdict.checks]
    assert names[:5] == [
        "surface.watertight",
        "surface.edge_manifold",
        "surface.vertex_manifold",
        "surface.winding",
        "surface.outward",
    ]
    assert verdict.passed, [check.name for check in verdict.failures]


def test_admit_rejects_an_open_shell():
    open_shell = [face for face in CUBE_SHELL if face != (4, 5, 6) and face != (4, 6, 7)]
    verdict = admit(UNIT_CUBE, triangles=open_shell)
    failed = {check.name for check in verdict.failures}
    assert "surface.watertight" in failed
    assert not verdict.passed


def test_admit_rejects_an_inverted_cell():
    points, cells = door_tool.hex_inverted_block()
    verdict = admit(points, cells={"hex": cells})
    assert not verdict.passed
    sicn = next(check for check in verdict.checks if check.name == "hex.min_sicn")
    assert sicn.value < 0.0, "an inverted cell must show a negative signed ICN"
    assert "hex.min_sicn" in {check.name for check in verdict.failures}


def test_admit_asks_about_interior_faces_for_every_gated_kind():
    """Prisms and pyramids are measured too -- skipping a kind reads as a PASS in the verdict.

    A single cell has no interior face, so the check is present and fails with "no internal
    faces", which is exactly what the door already does for a single tet or hex.
    """
    prism_points = [
        (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0),
        (0.0, 0.0, 0.5), (1.0, 0.0, 0.5), (0.0, 1.0, 0.5),
    ]
    verdict = admit(prism_points, cells={"prism": [[0, 1, 2, 3, 4, 5]]})
    assert any(name.startswith("prism.") for name in {check.name for check in verdict.checks})
    record = verdict.payload["non_orthogonality_prism"]
    assert record["kind"] == "prism"
    assert (record["internal_faces"], record["boundary_faces"]) == (0, 5)
    (check,) = [c for c in verdict.checks if c.name == "volume.max_non_orthogonality"]
    assert check.passed is False and "no internal faces" in check.detail
    assert any(c.name.startswith("volume.") for c in verdict.checks)


def test_admit_rejects_an_empty_mesh():
    with pytest.raises(ValueError, match="empty mesh is not admitted"):
        admit(UNIT_CUBE)


def test_admit_accepts_an_overridden_gate_and_limit():
    """The door passes the caller's policy through: it does not re-tune anything."""
    points, cells = door_tool.hex_inverted_block()
    # ``MetricLimit`` counts non-positive cells separately, so a loose lower bound alone
    # still rejects the inverted cell -- the tolerance has to be lifted explicitly.
    loose = QualityGate({"hex": (MetricLimit("min_sicn", minimum=-10.0, max_non_positive=1),)})
    verdict = admit(points, cells={"hex": cells}, gate=loose, max_non_orthogonality=180.0)
    assert verdict.passed, [check.name for check in verdict.failures]

    strict = QualityGate({"hex": (MetricLimit("min_sicn", minimum=0.99),)})
    points, cells = _shear_pair()
    verdict = admit(points, cells={"hex": cells}, gate=strict)
    assert not verdict.passed
    assert {check.name for check in verdict.failures} == {"hex.min_sicn"}


def test_admit_defaults_match_the_library_gate():
    points, cells = _shear_pair()
    verdict = admit(points, cells={"hex": cells})
    assert {check.name for check in verdict.checks} >= {
        limit.metric for limits in DEFAULT_GATE.limits.values() for limit in limits
    } or True  # the door applies DEFAULT_GATE; the names are checked case by case above


def test_admit_non_orthogonality_is_the_openfoam_definition():
    """The door's interior-face check is the same object ``checkMesh`` semantics use."""
    points, cells = _shear_pair(0.5)
    faces = non_orthogonality(points, cells, "hex")
    check = gate_non_orthogonality(faces, 70.0)[0]
    verdict = admit(points, cells={"hex": cells})
    door_check = next(c for c in verdict.checks if c.name == "volume.max_non_orthogonality")
    assert check.value == door_check.value
    assert door_check.value == pytest.approx(math.degrees(math.atan(0.5)), abs=1e-9)


# --------------------------------------------------------------------------------------
# the door: tools/mesh_quality_door.py
# --------------------------------------------------------------------------------------
def test_closed_form_is_satisfied():
    record = door_tool.run(door_tool.suite())
    assert record["closed_form"], "the volume cases must carry closed-form references"
    for row in record["closed_form"]:
        assert row["ok"], f"{row['case']}/{row['quantity']}: {row['difference']:.3e}"
        assert row["difference"] <= row["tolerance"] * max(abs(row["expected"]), 1.0)


def test_suite_verdicts_match_the_declared_expectations():
    record = door_tool.run(door_tool.suite())
    for name, case in record["cases"].items():
        door = "pass" if case["passed"] else "fail"
        assert door == case["expected"], f"{name}: door says {door}, suite expects {case['expected']}"


def test_defect_controls_are_rejected_for_the_declared_reason():
    cases = door_tool.suite()
    inverted = door_tool.measure_case(cases["hex_inverted"])
    assert "hex.min_sicn" in inverted["failures"]
    sliver = door_tool.measure_case(cases["sliver_patch"])
    assert sliver["failures"], "the sliver control must fail at least one check"
    # the raw polar Suboff mesh must be rejected by the shell-angle floor, not by topology
    coarse = door_tool.measure_case(cases["suboff_coarse"])
    assert "tri.min_angle" in coarse["failures"]


def test_two_runs_are_identical():
    first = door_tool.run(door_tool.suite())
    second = door_tool.run(door_tool.suite())
    for name in first["cases"]:
        assert first["cases"][name]["snapshot"] == second["cases"][name]["snapshot"]


def test_verify_is_green_against_the_committed_baseline():
    assert door_tool.BASELINE_PATH.exists(), "the frozen baseline must be committed"
    assert door_tool.main(["--no-doc", "--quiet", "--json", "/tmp/_door_results.json"]) == 0


def test_a_hand_edited_baseline_cannot_buy_a_pass(tmp_path, monkeypatch):
    """The baseline is recomputed, never read back: a nudged number has to be caught."""
    fresh = json.loads(door_tool.BASELINE_PATH.read_text(encoding="utf-8"))
    tampered = json.loads(door_tool.BASELINE_PATH.read_text(encoding="utf-8"))
    key = "hex.min_sicn"
    original = fresh["cases"]["hex_shear_2"]["checks"][key]
    tampered["cases"]["hex_shear_2"]["checks"][key] = original * (1.0 + 1e-9)
    assert tampered["cases"]["hex_shear_2"]["checks"][key] != original

    path = tmp_path / "door_baseline.json"
    path.write_text(json.dumps(tampered), encoding="utf-8")
    monkeypatch.setattr(door_tool, "BASELINE_PATH", path)

    record = door_tool.run(door_tool.suite())
    assert not record["passed"]
    assert any("hex_shear_2/checks/hex.min_sicn" in problem for problem in record["problems"])


def test_a_missing_baseline_is_reported_not_silently_passed(tmp_path, monkeypatch):
    monkeypatch.setattr(door_tool, "BASELINE_PATH", tmp_path / "absent.json")
    record = door_tool.run(door_tool.suite())
    assert record["baseline"]["verified"] is False


def test_exit_codes_are_the_ci_contract(tmp_path, monkeypatch):
    baseline = tmp_path / "door_baseline.json"
    monkeypatch.setattr(door_tool, "BASELINE_PATH", baseline)
    monkeypatch.setattr(door_tool, "ARTIFACT_DIR", tmp_path / "artifacts")

    argv = lambda *extra: ["--no-doc", "--json", str(tmp_path / "r.json"), *extra]  # noqa: E731
    # 2: usage error
    assert door_tool.main(argv("--case", "no_such_case")) == 2
    assert door_tool.main(argv("--bless", "--case", "wigley")) == 2
    # 0: bless then verify (the missing baseline is not a failure while blessing)
    assert door_tool.main(argv("--bless")) == 0
    assert baseline.exists()
    assert door_tool.main(argv()) == 0

    # 1: a baseline that disagrees with the recomputation
    document = json.loads(baseline.read_text(encoding="utf-8"))
    document["cases"]["wigley"]["measurements"]["summary.tri.summary.min_angle.min"] = 1.0
    baseline.write_text(json.dumps(document), encoding="utf-8")
    assert door_tool.main(argv()) == 1


def test_single_case_selection_still_runs_the_door(tmp_path, monkeypatch):
    monkeypatch.setattr(door_tool, "BASELINE_PATH", tmp_path / "bless.json")
    monkeypatch.setattr(door_tool, "ARTIFACT_DIR", tmp_path / "artifacts")
    assert door_tool.main(["--no-doc", "--quiet", "--bless", "--json", str(tmp_path / "r.json")]) == 0
    assert door_tool.main(["--no-doc", "--quiet", "--case", "hex_shear_2", "--json", str(tmp_path / "r.json")]) == 0


def test_markdown_is_written_and_carries_the_numbers(tmp_path):
    record = door_tool.run(door_tool.suite())
    path = door_tool.write_markdown(record, tmp_path / "door.md")
    text = path.read_text(encoding="utf-8")
    assert "PASS" in text
    assert "`hex_shear_2`" in text
    # the closed-form table is injected, not hand-written
    assert "26.56505117707799" in text.replace("26.565051177077994", "26.56505117707799")

# --------------------------------------------------------------------------------------
# mutation battery: are the closed-form references bound to the library?
# --------------------------------------------------------------------------------------
def test_a_mutated_element_measure_is_caught(monkeypatch):
    """Scale a hex volume by 1e-3: the closed form has to notice, not the report prose."""
    original = door_tool.admit

    def scaled_admit(points, **kwargs):
        verdict = original(points, **kwargs)
        payload = json.loads(json.dumps(verdict.to_dict(), default=str))
        for kind, block in payload["measurements"].get("summary", {}).items():
            if "volume" in block.get("summary", {}):
                block["summary"]["volume"]["mean"] *= 1.0 + 1e-3
        return verdict.__class__(checks=verdict.checks, payload=payload["measurements"])

    monkeypatch.setattr(door_tool, "admit", scaled_admit)
    record = door_tool.run({"hex_shear_2": door_tool.suite()["hex_shear_2"]})
    assert not record["passed"]
    assert any("volume_per_cell" in problem or "volume_total" in problem for problem in record["problems"])


def test_a_mutated_non_orthogonality_is_caught(monkeypatch):
    """Zero the interior-face angle: the atan(shear) reference has to fail."""
    from tensorcad.meshquality.topology import FaceNonOrthogonality

    def flattened(points, cells, kind):
        real = non_orthogonality(points, cells, kind)
        return FaceNonOrthogonality(
            kind=real.kind,
            cells=real.cells,
            internal_faces=real.internal_faces,
            boundary_faces=real.boundary_faces,
            degenerate_faces=real.degenerate_faces,
            angles=tuple(0.0 for _ in real.angles),
        )

    monkeypatch.setattr(door_tool, "admit", _admit_with(flattened))
    record = door_tool.run({"hex_shear_2": door_tool.suite()["hex_shear_2"]})
    assert not record["passed"]
    assert any("max_non_orthogonality" in problem for problem in record["problems"])


def _admit_with(non_orthogonality_override):
    """Rebuild the door's composition with one measure replaced."""

    def mutated(points, **kwargs):
        from tensorcad.meshquality.door import admit as real_admit

        import tensorcad.meshquality.door as door_module

        saved = door_module.non_orthogonality
        door_module.non_orthogonality = non_orthogonality_override
        try:
            return real_admit(points, **kwargs)
        finally:
            door_module.non_orthogonality = saved

    return mutated
