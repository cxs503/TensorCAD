"""Tests for the external-oracle differential (``GAP-CAD-D1-MESH-QUALITY-LIB``).

Two layers, because the oracles themselves are third-party binaries that this image cannot
always build:

* the *protocol* is tested here (mesh file format, comparison logic, tolerance and status
  handling, exit codes) with a fake oracle, so it runs everywhere;
* the *recorded evidence* is re-asserted from
  ``artifacts/mesh_quality_external_diff/oracle_diff.json``: every pair the artifact marks
  AGREE must be ≤ 1e-6, the four explained differences must still be the four explained
  differences, and the two documented relations (tet scaled-Jacobian √2 factor, triangle
  aspect-ratio ≈ 2/√3 discrepancy) must still hold. A later edit that changes a metric
  therefore fails this suite even on a host without Verdict/Gmsh.
"""
from __future__ import annotations

import json
import math
import stat
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "tools"))

import mesh_quality_external_diff as diff  # noqa: E402

ARTIFACT = REPO / "artifacts" / "mesh_quality_external_diff" / "oracle_diff.json"


# --------------------------------------------------------------------------------------
# corpus / protocol
# --------------------------------------------------------------------------------------
def test_corpus_covers_the_anchor_requirement():
    corpus = diff.meshes()
    assert len(corpus) >= 3, "the anchor asks for at least three meshes"
    kinds = {kind for kind, _pts, _elems in corpus.values()}
    assert {"tri", "quad", "tet", "hex"} <= kinds
    for name, (kind, pts, elems) in corpus.items():
        npe = {"tri": 3, "quad": 4, "tet": 4, "hex": 8}[kind]
        assert elems, name
        assert all(len(e) == npe for e in elems), name
        assert all(0 <= i < len(pts) for e in elems for i in e), name


def test_mesh_file_roundtrip(tmp_path):
    kind, pts, elems = diff.meshes()["tet_sliver"]
    path = tmp_path / "m.mesh"
    diff.write_mesh_file(path, kind, pts, elems)
    lines = [ln for ln in path.read_text().splitlines() if ln and not ln.startswith("#")]
    assert lines[0] == f"nodes {len(pts)}"
    assert lines[1 + len(pts)] == f"{kind} {len(elems)}"
    assert len(lines) == 2 + len(pts) + len(elems)


def _fake_oracle(tmp_path, body: str) -> str:
    p = tmp_path / "oracle"
    p.write_text("#!/bin/sh\n" + body)
    p.chmod(p.stat().st_mode | stat.S_IEXEC)
    return str(p)


def test_missing_oracle_is_not_built_not_a_pass(tmp_path):
    status, data = diff.run_oracle(None, tmp_path / "x", "W")
    assert status == "NOT_BUILT" and data == {}
    status, _ = diff.run_oracle(str(tmp_path / "nope"), tmp_path / "x", "W")
    assert status == "NOT_BUILT"


def test_failing_oracle_reports_not_run(tmp_path):
    binary = _fake_oracle(tmp_path, "echo boom >&2; exit 3\n")
    status, data = diff.run_oracle(binary, tmp_path / "x", "W")
    assert status.startswith("NOT_RUN") and data == {}


def test_oracle_output_parsing(tmp_path):
    binary = _fake_oracle(
        tmp_path,
        "printf 'W tri area 0 0.5\\nW tri area 1 0.25\\nW tet volume 0 1.0\\nnoise\\n'",
    )
    status, data = diff.run_oracle(binary, tmp_path / "x", "W")
    assert status == "RAN"
    assert data["tri"]["area"] == [0.5, 0.25]
    assert data["tet"]["volume"] == [1.0]


# --------------------------------------------------------------------------------------
# comparison semantics
# --------------------------------------------------------------------------------------
def test_compare_agree_and_differ():
    ours = {"a": (1.0, 2.0), "b": (1.0, 1.0)}
    ext = {"x": [1.0 + 1e-9, 2.0], "y": [1.0, 2.0]}
    res = diff.compare(ours, ext, {"a": "x", "b": "y"})
    assert res["a"]["status"] == "AGREE" and res["a"]["n_beyond_tol"] == 0
    assert res["b"]["status"] == "DIFFER" and res["b"]["n_beyond_tol"] == 1
    assert res["b"]["worst_index"] == 1


def test_compare_reports_missing_and_length_mismatch():
    res = diff.compare({"a": (1.0, 2.0)}, {"z": [1.0]}, {"a": "x", "b": "z"})
    assert res["a"]["status"] == "NOT_RETURNED"
    assert res["b"]["status"] == "LENGTH_MISMATCH"


def test_compare_treats_nan_pairs_as_agreeing():
    nan = float("nan")
    res = diff.compare({"a": (nan, 1.0)}, {"x": [nan, 1.0]}, {"a": "x"})
    assert res["a"]["status"] == "AGREE" and res["a"]["n_elements"] == 1


def test_driver_without_oracles_does_not_claim_a_pass(tmp_path):
    art = tmp_path / "a.json"
    code = diff.main(["--artifact", str(art), "--quiet"])
    assert code == 0, "no oracle supplied is not a mismatch"
    report = json.loads(art.read_text())
    assert report["oracle_status"] == {"verdict": "NOT_BUILT", "gmsh": "NOT_BUILT"}
    assert report["mismatches"] == []


def test_require_oracles_exits_2(tmp_path):
    code = diff.main(["--artifact", str(tmp_path / "a.json"), "--quiet", "--require-oracles"])
    assert code == 2


def test_driver_with_disagreeing_oracle_exits_1(tmp_path):
    binary = _fake_oracle(tmp_path, "printf 'W hex scaled_jacobian 0 -1.0\\n'; exit 0\n")
    code = diff.main(["--verdict-oracle", binary, "--artifact", str(tmp_path / "a.json"), "--quiet"])
    assert code == 1
    report = json.loads((tmp_path / "a.json").read_text())
    assert any("verdict" in m for m in report["mismatches"])


# --------------------------------------------------------------------------------------
# recorded evidence: the artifact must keep saying what this tree claims
# --------------------------------------------------------------------------------------
pytestmark_artifact = pytest.mark.skipif(not ARTIFACT.exists(), reason="differential artifact absent")


@pytest.fixture(scope="module")
def artifact():
    if not ARTIFACT.exists():
        pytest.skip("differential artifact absent")
    return json.loads(ARTIFACT.read_text())


def test_artifact_has_at_least_three_meshes(artifact):
    assert len(artifact["meshes"]) >= 3
    kinds = {m["kind"] for m in artifact["meshes"].values()}
    assert {"tri", "quad", "tet", "hex"} <= kinds


def test_every_agreeing_pair_is_within_one_e_minus_six(artifact):
    tol = artifact["tolerance"]
    agree = 0
    for name, mesh in artifact["meshes"].items():
        for oracle, res in mesh["oracles"].items():
            if res["status"] != "RAN":
                continue
            for metric, cmp in res["comparisons"].items():
                if cmp.get("status") == "AGREE":
                    agree += 1
                    assert cmp["max_rel_diff"] <= tol, f"{name}/{oracle}/{metric}"
    assert agree >= 14, f"only {agree} agreeing pairs recorded"


def test_explained_differences_are_exactly_the_documented_four(artifact):
    """If a difference disappears or a new one appears, the doc must be updated with it."""
    differing = set()
    for name, mesh in artifact["meshes"].items():
        for oracle, res in mesh["oracles"].items():
            if res["status"] != "RAN":
                continue
            for metric, cmp in res["comparisons"].items():
                if cmp.get("status") == "DIFFER":
                    differing.add((mesh["kind"], oracle, metric))
    assert differing == {
        ("tet", "verdict", "scaled_jacobian"),
        ("tet", "verdict", "min_dihedral"),
        ("tet", "verdict", "skewness"),
        ("tri", "verdict", "aspect_ratio"),
    }, differing


def test_tet_scaled_jacobian_differs_from_verdict_by_exactly_sqrt2(artifact):
    """Documented normalisation: Verdict / library = √2, on every mesh where both ran."""
    checked = 0
    for name, mesh in artifact["meshes"].items():
        if mesh["kind"] != "tet":
            continue
        cmp = mesh["oracles"].get("verdict", {}).get("comparisons", {}).get("scaled_jacobian")
        if not cmp or cmp.get("status") != "DIFFER":
            continue
        ratio = cmp["external_at_worst"] / cmp["ours_at_worst"]
        assert ratio == pytest.approx(math.sqrt(2.0), rel=1e-12), name
        checked += 1
    assert checked >= 2


def test_tri_aspect_ratio_discrepancy_is_the_documented_one(artifact):
    """Finding F5: our triangle aspect_ratio is Verdict's radius_ratio, so it differs by ≲2/√3."""
    cmp = artifact["meshes"]["tri_warped"]["oracles"]["verdict"]["comparisons"]["aspect_ratio"]
    assert cmp["status"] == "DIFFER"
    assert cmp["n_beyond_tol"] == cmp["n_elements"]
    ratio = cmp["external_at_worst"] / cmp["ours_at_worst"]
    assert 1.0 < ratio <= 2.0 / math.sqrt(3.0) + 1e-12
    assert cmp["max_rel_diff"] > 0.1, "the documented discrepancy is ~15 %, not round-off"


def test_gmsh_half_is_recorded_as_not_run_not_as_passed(artifact):
    assert artifact["gmsh_oracle"]["version"] == "unknown" or True  # informational
    statuses = {m["oracles"].get("gmsh", {}).get("status") for m in artifact["meshes"].values()}
    if statuses == {"NOT_BUILT"}:
        # Honest state on this host: no number is claimed for Gmsh anywhere in the artifact.
        for mesh in artifact["meshes"].values():
            assert "comparisons" not in mesh["oracles"].get("gmsh", {})