"""Single-body 10^6-10^7 driver: artifact shape, caps and drift binding.

Two failure modes are what these tests defend against, both of which have bitten
scale benchmarks in this repository before:

* a rung that was *refused* (too large) silently disappearing, so the table looks
  complete while the biggest size was never measured -- ``not_run`` must be
  recorded with a reason;
* the document drifting away from the committed measurement JSON -- the doc block
  must be exactly what the artifact renders to.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from tensorcad.perf_evidence import (
    SINGLE10M_TABLE_NAMES,
    check_document,
    extract_block,
    render_block,
    render_document,
)
from tensorcad.perf_subdiv import run_single_body_benchmark

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "benchmarks" / "perf_scale" / "artifacts" / "ladder_single_10m.json"
DOC = ROOT / "docs" / "perf_scale_single_body_10m.md"


@pytest.fixture(scope="module")
def small_artifact():
    return run_single_body_benchmark(
        bases=[(41, 13)],
        levels_ladder=[0, 1, 2],
        n_queries=32,
        max_triangles=10_000_000,
    )


def test_rows_cover_every_level_and_subtask(small_artifact):
    results = small_artifact["results"]
    subtasks = {r["subtask"] for r in results}
    assert subtasks == {"load", "subdivide", "volume", "point_locate_accel"}
    levels = {r["levels"] for r in results if r["subtask"] == "volume"}
    assert levels == {0, 1, 2}
    # level 0 has no subdivision step: the hull *is* the rung.
    assert [r["levels"] for r in results if r["subtask"] == "subdivide"] == [1, 2]
    assert small_artifact["schema"] == "tensorcad.perf_single10m/1"


def test_invariants_recorded_and_passing(small_artifact):
    summary = small_artifact["summary"]
    assert summary["all_watertight"] and summary["all_edge_manifold"]
    assert summary["all_euler_2"]
    assert summary["max_volume_rel_drift_vs_base"] < 1e-12
    assert summary["max_surface_area_rel_drift"] < 1e-12
    assert summary["invariants_ok"] is True
    # The analytic error is a parent-discretisation property, constant on every
    # rung -- it must not be folded into the invariant verdict.
    assert summary["max_analytic_rel_err"] == pytest.approx(0.00880431, abs=1e-6)


def test_volume_matches_analytic_reference(small_artifact):
    exact = small_artifact["hull"]["analytic_volume_m3"]
    for c in small_artifact["checks"]:
        assert c["analytic_volume_m3"] == pytest.approx(exact)
        assert abs(c["volume_rel_err"]) < 0.02


def test_ray_counts_and_index_are_reported(small_artifact):
    rows = [r for r in small_artifact["results"] if r["subtask"] == "point_locate_accel"]
    assert len(rows) == 3
    for row in rows:
        assert row["status"] == "measured"
        assert row["seconds"] > 0
        assert row["index_cells"] > 0 and row["index_bytes"] > 0
        assert "build=" in row["note"]


def test_size_cap_records_not_run_instead_of_dropping():
    artifact = run_single_body_benchmark(
        bases=[(41, 13)], levels_ladder=[0, 1, 2, 3], n_queries=16,
        max_triangles=8_000,
    )
    refused = [r for r in artifact["results"] if r["status"] == "not_run"]
    assert refused, "a rung above the cap must be recorded, not skipped"
    assert all("max_triangles" in r["note"] for r in refused)
    assert artifact["summary"]["single_body_10m_reached"] is False


def test_renderers_are_stable(small_artifact):
    artifacts = {"single10m": small_artifact}
    for name in SINGLE10M_TABLE_NAMES:
        body = render_block(name, artifacts)
        assert body.strip()
        assert body.count("|") > 4
    table = render_block("single10m-table", artifacts)
    assert "31,968" in table  # level-2 rung present
    exponents = render_block("single10m-exponents", artifacts)
    assert "not_run" not in exponents or "规模点不足" in exponents


def test_document_drift_is_detected(small_artifact):
    artifacts = {"single10m": small_artifact}
    body = render_block("single10m-table", artifacts)
    doc = (
        "text\n<!-- perf-evidence:single10m-table begin -->\n"
        + body
        + "\n<!-- perf-evidence:single10m-table end -->\n"
    )
    assert check_document(doc, artifacts, names=("single10m-table",)) == []
    tampered = doc.replace("31,968", "31,969")
    problems = check_document(tampered, artifacts, names=("single10m-table",))
    assert problems and "differs" in problems[0]


def test_document_missing_marker_is_drift(small_artifact):
    artifacts = {"single10m": small_artifact}
    problems = check_document("no tables here", artifacts, names=("single10m-table",))
    assert problems and "markers not found" in problems[0]


def test_render_document_rewrites_block(small_artifact):
    artifacts = {"single10m": small_artifact}
    doc = (
        "<!-- perf-evidence:single10m-table begin -->\nstale\n"
        "<!-- perf-evidence:single10m-table end -->\n"
    )
    out = render_document(doc, artifacts, names=("single10m-table",))
    assert "stale" not in out
    assert extract_block(out, "single10m-table") == render_block("single10m-table", artifacts)


@pytest.mark.skipif(not ARTIFACT.exists(), reason="10^6-10^7 single-body artifact not committed")
def test_committed_document_matches_committed_artifact():
    """The published 10^6-10^7 table must be rendered from the committed JSON."""
    from tensorcad.perf_evidence import load_artifact

    artifacts = {"single10m": load_artifact(ARTIFACT)}
    text = DOC.read_text(encoding="utf-8")
    problems = check_document(text, artifacts, names=SINGLE10M_TABLE_NAMES)
    assert problems == []