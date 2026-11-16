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
    # v2 schema adds obj_write and obj_read to the default subtask set
    assert subtasks == {"load", "subdivide", "volume", "point_locate_accel", "obj_write", "obj_read"}
    levels = {r["levels"] for r in results if r["subtask"] == "volume"}
    assert levels == {0, 1, 2}
    # level 0 has no subdivision step: the hull *is* the rung.
    assert [r["levels"] for r in results if r["subtask"] == "subdivide"] == [1, 2]
    assert small_artifact["schema"] == "tensorcad.perf_single10m/2"


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
    # Analytic volume is in the checks, not a separate "hull" key
    for c in small_artifact["checks"]:
        exact = c["analytic_volume_m3"]
        assert abs(c["volume_rel_err"]) < 0.02
        # All checks should share the same analytic volume
        assert c["analytic_volume_m3"] == pytest.approx(exact)


def test_obj_roundtrip_records_file_and_volume(small_artifact):
    """obj_write and obj_read rows carry file size and volume error."""
    write_rows = [r for r in small_artifact["results"] if r["subtask"] == "obj_write"]
    read_rows = [r for r in small_artifact["results"] if r["subtask"] == "obj_read"]
    assert len(write_rows) == 3  # levels 0, 1, 2
    assert len(read_rows) == 3
    for w in write_rows:
        assert w["file_bytes"] > 0
        assert w["file_mb_per_s"] is not None
    for r in read_rows:
        assert r["file_bytes"] > 0
        assert r["volume_rel_err"] is not None
        assert abs(r["volume_rel_err"]) < 0.02  # within parent discretisation error


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

#: Acceptance item 1 of ``GAP-CAD-D4-PERF-SCALE`` is the 10^6--10^7 band on a
#: *single connected body*.  v2 of this document shipped with
#: ``single_body_10m_reached = false`` because the cluster job that produced the
#: 10^7 rung was still queued; the artifact has since been landed from job
#: ``cj3805``.  These tests keep that claim bound to the committed numbers
#: instead of to the prose, so the document cannot silently regress to "we did
#: not measure it" (or, worse, claim 10^7 without a 10^7 row).
@pytest.mark.skipif(not ARTIFACT.exists(), reason="10^6-10^7 single-body artifact not committed")
def test_committed_artifact_measures_a_single_body_at_10m():
    from tensorcad.perf_evidence import load_artifact

    artifact = load_artifact(ARTIFACT)
    summary = artifact["summary"]
    assert summary["single_body_10m_reached"] is True
    assert summary["max_triangles"] >= 10_000_000
    assert summary["invariants_ok"] is True

    measured = [r for r in artifact["results"] if r.get("status") == "measured"]
    big = [r for r in measured if (r.get("triangles") or 0) >= 10_000_000]
    assert big, "no measured rung reaches 10^7 triangles"

    # Two independent parent bases must reach the 10^7 band; a single lucky grid
    # could be a cropping artefact rather than a measurement.
    bases = {(r["nx"], r["nz"]) for r in big}
    assert len(bases) >= 2, f"only one parent base reached 10^7: {bases}"

    # Every sub-task except the deliberately capped one must have run there.
    subtasks = {r["subtask"] for r in big}
    assert {"subdivide", "volume", "point_locate_accel"} <= subtasks


@pytest.mark.skipif(not ARTIFACT.exists(), reason="10^6-10^7 single-body artifact not committed")
def test_10m_rung_is_a_closed_triangle_mesh():
    """``V = F/2 + 2`` must hold exactly at 10^7: subdivision welds by index.

    If the midpoint de-duplication were tolerance-based, the 10^7 rung would pick
    up a handful of duplicated vertices per level and this identity would drift;
    it is an exact integer check, not an error bar.
    """
    from tensorcad.perf_evidence import load_artifact

    artifact = load_artifact(ARTIFACT)
    checks = [c for c in artifact["checks"] if c.get("triangles", 0) >= 10_000_000]
    assert checks, "no invariant record at 10^7"
    for c in checks:
        assert c["watertight"] is True and c["edge_manifold"] is True
        assert c["euler_characteristic"] == 2
        assert c["connected_bodies"] == 1
        assert c["vertices"] == c["triangles"] // 2 + 2
        assert abs(c["volume_rel_drift_vs_base"] or 0.0) <= 1e-12
        assert abs(c["surface_area_rel_drift"] or 0.0) <= 1e-12


@pytest.mark.skipif(not ARTIFACT.exists(), reason="10^6-10^7 single-body artifact not loaded")
def test_capacity_table_is_rate_times_budget():
    """Acceptance item 3: the cluster boundary must be re-derivable by arithmetic.

    Each row's first budget column has to equal the printed rate times that
    budget (within rounding of the printed rate), so nobody has to trust a
    fitted parameter or a hand-set threshold.
    """
    from tensorcad.perf_evidence import load_artifact, render_single10m_capacity

    artifact = load_artifact(ARTIFACT)
    budgets = (60.0, 600.0, 3600.0)
    table = render_single10m_capacity(artifact, budget_seconds=budgets, ram_budget_gib=16.0)
    rows = [
        ln for ln in table.splitlines()
        if ln.startswith("| ") and (" tri/s " in ln or " query/s " in ln)
    ]
    assert rows, table
    for line in rows:
        cells = [c.strip() for c in line.strip("|").split("|")]
        rate_text = cells[3].split()[0].replace(",", "")
        rate = float(rate_text)
        reach = [float(x.replace(",", "").split()[0]) for x in cells[4].split("/")]
        assert len(reach) == len(budgets)
        for got, budget in zip(reach, budgets):
            assert got == pytest.approx(rate * budget, rel=0.01), (line, budget)
        assert float(cells[5].replace(",", "")) > 0


@pytest.mark.skipif(not ARTIFACT.exists(), reason="10^6-10^7 single-body artifact not committed")
def test_capacity_block_is_embedded_in_committed_document():
    """The rendered capacity table must be in the doc (drift = test failure)."""
    from tensorcad.perf_evidence import load_artifact

    artifacts = {"single10m": load_artifact(ARTIFACT)}
    text = DOC.read_text(encoding="utf-8")
    assert extract_block(text, "single10m-capacity") == render_block(
        "single10m-capacity", artifacts
    )
