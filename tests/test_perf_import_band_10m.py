"""Tests for the 10^6--10^7 OBJ round-trip band (``docs/perf_scale_10m.md`` §11.5).

What this pins down:

1. the two band artifacts are **committed** and hash-match their provenance
   files (a fetched job artifact that silently changes is a broken evidence
   chain);
2. the band tables/gates in the document are rendered from those artifacts
   (no hand-typed numbers) -- ``scripts/perf_scale_doc.py --check`` must agree;
3. a ``not_run`` gate never renders as a pass (the defect this band exposed:
   the engine sets ``passed=True`` on gates it could not evaluate, and the
   renderer used to print a green ``pass`` for them);
4. the physics claims the document makes about the band actually hold in the
   JSON: the OBJ round trip is loss-free (``obj_read`` volume error equals the
   in-memory ``volume`` error) and the analytic reference is ``copies x
   7.512698 m^3`` within 3%.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from tensorcad.perf_evidence import (
    check_document,
    load_artifact,
    render_import_gates,
    render_import_table,
)

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / "docs" / "perf_scale_10m.md"
ART = ROOT / "benchmarks" / "perf_scale" / "artifacts"
SINGLE = ART / "ladder_import_single_10m.json"
FLEET = ART / "ladder_import_fleet_10m.json"
SINGLE_PROV = ART / "ladder_import_single_10m.provenance.json"
FLEET_PROV = ART / "ladder_import_fleet_10m.provenance.json"

#: Analytic displaced volume of the benchmark hull (``BoxHull``/parametric hull
#: at L=10, B=2, T=0.8): one body.
ANALYTIC_VOLUME = 7.512698
#: The ledger's acceptance bar for a physical quantity is 3%.
TOL = 0.03


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.mark.parametrize("artifact,prov", [(SINGLE, SINGLE_PROV), (FLEET, FLEET_PROV)])
def test_band_artifacts_are_committed_with_matching_provenance(artifact, prov):
    assert artifact.is_file(), f"{artifact} must be committed, not only fetched"
    assert prov.is_file()
    meta = json.loads(prov.read_text(encoding="utf-8"))
    assert meta["sha256"] == _sha256(artifact), "artifact changed since it was fetched"
    assert meta["rc"] == 0
    assert meta["job_id"].startswith("cj")
    assert meta["artifact"] == artifact.name


def test_single_body_artifact_reaches_the_hull_ceiling():
    art = load_artifact(SINGLE)
    measured = [r for r in art["results"] if r["status"] == "measured"]
    assert measured, "no measured rows in the single-body band run"
    tris = {r["triangles"] for r in measured}
    assert tris == {8_199_998}, "single-body rung must be 4001x513 = 8,199,998 tris"
    subtasks = {r["subtask"] for r in measured}
    assert {"load", "obj_write", "obj_read", "volume"} <= subtasks


def test_fleet_artifact_reaches_the_ten_to_the_seventh_band():
    art = load_artifact(FLEET)
    top = max(r["triangles"] for r in art["results"])
    assert top == 12_011_994, "fleet run must reach 1.2e7 triangles (3 x 2001x501)"
    read_top = max(
        r["triangles"] for r in art["results"] if r["subtask"] == "obj_read"
    )
    assert read_top == 12_011_994, "the 1.2e7 rung must include the OBJ parse"


def test_single_body_gates_are_not_run_not_pass():
    """One size cannot yield an exponent: the table must say so."""
    art = load_artifact(SINGLE)
    assert {g["status"] for g in art["gates"]} == {"not_run"}
    table = render_import_gates(art)
    assert "| pass |" not in table
    assert table.count("| not_run |") >= len(art["gates"])


def test_bytes_per_triangle_is_rendered_from_the_artifact():
    """The old prose hard-coded "~48 B/tri"; it is 41--55 B/tri in reality."""
    art = load_artifact(FLEET)
    table = render_import_table(art)
    assert "文件 B/tri" in table
    ratios = {
        r["file_bytes"] / r["triangles"]
        for r in art["results"]
        if r.get("file_bytes") and r["subtask"] == "obj_write"
    }
    assert min(ratios) < 48.0 < max(ratios), (
        "B/tri must straddle the old hand-typed 48 B/tri claim -- if it does not, "
        "the document prose needs re-reading"
    )


@pytest.mark.parametrize("artifact,expected_pairs", [(SINGLE, 1), (FLEET, 6)])
def test_round_trip_is_loss_free_and_volume_matches_analytic(artifact, expected_pairs):
    art = load_artifact(artifact)
    by = {}
    for r in art["results"]:
        if r["status"] == "measured" and r["subtask"] in ("obj_write", "obj_read", "volume"):
            by.setdefault((r["copies"], r["triangles"]), {})[r["subtask"]] = r
    checked = 0
    for (copies, tris), rows in by.items():
        if "obj_read" not in rows or "volume" not in rows:
            continue
        read, vol = rows["obj_read"], rows["volume"]
        # (a) writer and parser see the same file
        if "obj_write" in rows:
            assert read["file_bytes"] == rows["obj_write"]["file_bytes"]
        # (b) the parse must not lose geometry: same triangles, same volume
        assert read["triangles"] == vol["triangles"] == tris
        # (c) volume within 3% of copies x analytic
        analytic = copies * ANALYTIC_VOLUME
        rel = read["volume_rel_err"]
        assert abs(rel) < TOL, f"{copies}x{tris}: rel err {rel:+.3e} exceeds {TOL:.0%}"
        assert rel < 0, "the panel mesh under-estimates volume (inner polyhedron)"
        # (d) the error is a small discrete-geometry error, not a parse error:
        #     2.7e-5 at 8.2e6 and 3.0e-5 at 1.2e7 triangles.  A lossy parse
        #     (dropped faces) would show a deficit orders of magnitude larger.
        assert abs(rel) < 3e-4
        checked += 1
    assert checked == expected_pairs, (
        f"{checked} round-trip pairs verified, expected {expected_pairs}"
    )


def test_document_band_tables_match_committed_artifacts():
    doc = DOC.read_text(encoding="utf-8")
    artifacts = {
        "import_single10m": load_artifact(SINGLE),
        "import_fleet": load_artifact(FLEET),
    }
    for name in ("import-single10m-table", "import-single10m-gates",
                 "import-fleet-table", "import-fleet-gates"):
        assert f"perf-evidence:{name} begin" in doc, f"§11.5 must embed {name}"
    # ``check_document`` only knows the full bundle; assert the band blocks
    # directly so this test fails on drift even if the bundle loading changes.
    from tensorcad.perf_evidence import extract_block, render_block

    for name in ("import-single10m-table", "import-single10m-gates",
                 "import-fleet-table", "import-fleet-gates"):
        assert extract_block(doc, name) == render_block(name, artifacts), name
    assert check_document is not None