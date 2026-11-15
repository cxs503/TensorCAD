"""Tests for the large-geometry evidence bundle (``tensorcad.perf_evidence``).

Two things are asserted here:

1. the renderers are pure functions of the measurement JSON (no hand-typed
   numbers, ``not_run`` rows stay explicitly ``not_run``);
2. the *committed* document and the *committed* artifacts agree -- so the
   published scaling table cannot drift from the run that produced it.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from tensorcad.perf_evidence import (
    DASH,
    TABLE_NAMES,
    check_document,
    extract_block,
    load_artifact,
    marker_begin,
    marker_end,
    render_block,
    render_document,
    render_budget_table,
    render_fleet3pt_gates,
    render_fleet3pt_table,
    render_mesh_table,
    render_voxel_table,
)

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / "docs" / "perf_scale_10m.md"
MESH = ROOT / "benchmarks" / "perf_scale" / "artifacts" / "ladder_mesh.json"
VOXEL = ROOT / "benchmarks" / "perf_scale" / "artifacts" / "ladder_voxel.json"
FLEET3PT = ROOT / "benchmarks" / "perf_scale" / "ladder_fleet_10m_3pt.json"


def _row(**kw):
    base = {
        "subtask": "load", "copies": 1, "nx": 161, "nz": 41,
        "triangles": 25918, "vertices": 12961, "grid_cells": 0,
        "seconds": 1.2345, "status": "measured", "throughput": 20998.0,
        "note": "", "build_seconds": None, "index_cells": 0,
        "index_entries": 0, "index_bytes": 0,
        "peak_rss_mb": 640.4, "rss_mb": 120.0, "mesh_mb": 9.5,
    }
    base.update(kw)
    return base


def _mesh_artifact(rows):
    return {"schema": "tensorcad.perf_benchmark/2", "results": rows}


def _fleet3pt_artifact(rows=None, gates=None, checks=None):
    return {
        "schema": "tensorcad.perf_benchmark/2",
        "results": rows if rows is not None else [_row()],
        "gates": gates or [],
        "gates_checked": len(gates or []),
        "gates_not_run": sum(1 for g in (gates or []) if g["status"] == "not_run"),
        "volume_checks": checks or [],
        "verdict": "pass" if gates and all(g["passed"] for g in gates) else "not_run",
    }


def _all(mesh, voxel, fleet3pt=None):
    """Artifact mapping that satisfies every marked table in the document."""
    return {"mesh": mesh, "voxel": voxel,
            "fleet3pt": fleet3pt if fleet3pt is not None else _fleet3pt_artifact()}


# ---------------------------------------------------------------------------
# Renderers
# ---------------------------------------------------------------------------

def test_mesh_table_has_one_row_per_measurement() -> None:
    rows = [
        _row(subtask="volume", triangles=4003998, seconds=1.644, throughput=2_435_000.0),
        _row(subtask="load", triangles=678, seconds=0.019, throughput=34962.0),
        _row(subtask="point_locate_accel", triangles=4003998, seconds=1.09,
             throughput=1836.0, build_seconds=31.11, index_bytes=334_900_000,
             peak_rss_mb=6151.0),
    ]
    table = render_mesh_table(_mesh_artifact(rows))
    lines = table.splitlines()
    assert lines[0].startswith("| 子任务")
    # header + separator + one line per result
    assert len(lines) == 2 + len(rows)
    # sub-task order follows the display order (load first, volume last)
    assert "load" in lines[2]
    assert "point_locate_accel" in lines[3]
    assert "volume" in lines[4]


def test_mesh_table_formats_index_cost_and_memory() -> None:
    table = render_mesh_table(_mesh_artifact([
        _row(subtask="point_locate_accel", build_seconds=31.114, index_bytes=334_900_000),
    ]))
    body = table.splitlines()[2]
    assert "31.11" in body            # build seconds
    assert "319.4" in body            # 334,900,000 B -> MiB
    assert "640" in body              # peak RSS


def test_mesh_table_keeps_not_run_explicit() -> None:
    table = render_mesh_table(_mesh_artifact([
        _row(subtask="point_locate", triangles=1_640_958, seconds=None,
             throughput=None, status="not_run", note="超阈值跳过"),
    ]))
    body = table.splitlines()[2]
    assert "not_run" in body
    assert "超阈值跳过" in body
    assert DASH in body


def test_voxel_table_sorted_by_cells() -> None:
    rows = [
        _row(subtask="voxelize", nx=384, nz=0, triangles=0, vertices=0,
             grid_cells=56_623_104, seconds=366.12, throughput=154658.0),
        _row(subtask="voxelize", nx=32, nz=0, triangles=0, vertices=0,
             grid_cells=32768, seconds=0.287, throughput=114029.0),
    ]
    lines = render_voxel_table(_mesh_artifact(rows)).splitlines()
    assert "32³" in lines[2]
    assert "384³" in lines[3]


def test_render_block_rejects_unknown_table() -> None:
    with pytest.raises(ValueError):
        render_block("no-such-table", _mesh_artifact([_row()]))


def test_render_is_deterministic() -> None:
    artifact = _mesh_artifact([_row(), _row(subtask="volume")])
    assert render_mesh_table(artifact) == render_mesh_table(artifact)


# ---------------------------------------------------------------------------
# Document binding
# ---------------------------------------------------------------------------

def _doc_with(blocks: dict[str, str]) -> str:
    parts = ["# doc", ""]
    for name, body in blocks.items():
        parts.append(f"## {name}")
        parts.append(marker_begin(name))
        parts.append(body)
        parts.append(marker_end(name))
    return "\n".join(parts) + "\n"


def test_check_document_accepts_matching_doc() -> None:
    mesh = _mesh_artifact([_row()])
    voxel = _mesh_artifact([_row(subtask="voxelize", nx=32, nz=0, triangles=0,
                                grid_cells=32768)])
    artifacts = _all(mesh, voxel)
    doc = _doc_with({name: render_block(name, artifacts) for name in TABLE_NAMES})
    assert check_document(doc, _all(mesh, voxel)) == []


def test_check_document_detects_hand_edited_number() -> None:
    mesh = _mesh_artifact([_row(seconds=1.2345)])
    voxel = _mesh_artifact([_row(subtask="voxelize", nx=32, nz=0, triangles=0,
                                grid_cells=32768)])
    artifacts = _all(mesh, voxel)
    blocks = {name: render_block(name, artifacts) for name in TABLE_NAMES}
    blocks["mesh-table"] = blocks["mesh-table"].replace("1.234", "0.987")
    doc = _doc_with(blocks)
    problems = check_document(doc, _all(mesh, voxel))
    assert len(problems) == 1
    assert "mesh-table" in problems[0]


def test_check_document_reports_missing_markers() -> None:
    mesh = _mesh_artifact([_row()])
    problems = check_document("# no tables here\n", {"mesh": mesh})
    assert any("markers not found" in p for p in problems)


def test_render_document_is_idempotent() -> None:
    mesh = _mesh_artifact([_row()])
    voxel = _mesh_artifact([_row(subtask="voxelize", nx=32, nz=0, triangles=0,
                                grid_cells=32768)])
    artifacts = _all(mesh, voxel)
    once = render_document(_doc_with({k: "" for k in TABLE_NAMES}), artifacts)
    twice = render_document(once, artifacts)
    assert once == twice
    assert check_document(once, artifacts) == []


def test_budget_table_extrapolates_measured_rate() -> None:
    mesh = _mesh_artifact([
        _row(subtask="load", triangles=4003998, seconds=213.7, throughput=18737.0),
    ])
    voxel = _mesh_artifact([
        _row(subtask="voxelize", nx=384, nz=0, triangles=0, vertices=0,
             grid_cells=56_623_104, seconds=366.12, throughput=154658.0),
    ])
    lines = render_budget_table(mesh, voxel).splitlines()
    assert len(lines) == 2 + 2
    load_row = [ln for ln in lines if "| load |" in ln][0]
    # 18,737 tri/s * 60 s
    assert "1,124,220" in load_row
    assert "213.700" in load_row
    voxel_row = [ln for ln in lines if "| voxelize |" in ln][0]
    assert "9,279,480" in voxel_row       # 154,658 cell/s * 60 s


# ---------------------------------------------------------------------------
# Committed evidence (integration)
# ---------------------------------------------------------------------------

def test_committed_artifacts_exist_and_are_measurements() -> None:
    for path in (MESH, VOXEL):
        assert path.exists(), f"missing committed artifact {path}"
        artifact = load_artifact(path)
        rows = artifact["results"]
        assert rows, f"{path} has no results"
        assert any(r["status"] == "measured" for r in rows)
        assert artifact.get("environment"), f"{path} has no environment record"


def test_committed_document_matches_committed_artifacts() -> None:
    assert DOC.exists(), f"missing {DOC}"
    problems = check_document(
        DOC.read_text(encoding="utf-8"),
        {"mesh": load_artifact(MESH), "voxel": load_artifact(VOXEL),
         "fleet3pt": load_artifact(FLEET3PT)},
    )
    assert problems == [], (
        "docs/perf_scale_10m.md drifted from the committed artifacts; "
        "rerun `python scripts/perf_scale_doc.py --write`"
    )


def test_mesh_artifact_reaches_the_1e6_band() -> None:
    """The anchor is about 10^6--10^7 triangles: the evidence must get there."""
    mesh = load_artifact(MESH)
    biggest = max(r["triangles"] for r in mesh["results"])
    assert biggest >= 1_000_000, f"largest rung is only {biggest:,} triangles"


def test_accelerated_point_location_beats_brute_force_at_same_size() -> None:
    """The published speed-up claim, re-derived from the artifact."""
    mesh = load_artifact(MESH)
    by_size: dict[int, dict[str, float]] = {}
    for r in mesh["results"]:
        if r["status"] != "measured" or not r["seconds"]:
            continue
        by_size.setdefault(r["triangles"], {})[r["subtask"]] = r["seconds"]
    compared = 0
    for size, times in by_size.items():
        if "point_locate" in times and "point_locate_accel" in times:
            assert times["point_locate_accel"] < times["point_locate"], size
            compared += 1
    assert compared >= 3, "not enough sizes with both point-location paths"


# ---------------------------------------------------------------------------
# 10^6-10^7 multi-body ladder (fleet3pt)
# ---------------------------------------------------------------------------

def _gate(name, detail, passed=True, status="measured", ratio=None):
    return {"name": name, "detail": detail, "passed": passed,
            "status": status, "ratio": ratio}


def test_fleet3pt_table_one_row_per_rung_in_subtask_order() -> None:
    rows = [
        _row(subtask="volume", copies=3, triangles=12_011_994),
        _row(subtask="load", copies=1, triangles=4_003_998),
        _row(subtask="point_locate_accel", copies=2, triangles=8_007_996,
             build_seconds=63.55, index_bytes=636_834_220),
        _row(subtask="fleet_assemble", copies=2, triangles=8_007_996),
    ]
    lines = render_fleet3pt_table(_fleet3pt_artifact(rows)).splitlines()
    assert len(lines) == 2 + len(rows)
    assert "load" in lines[2] and "fleet_assemble" in lines[3]
    assert "point_locate_accel" in lines[4] and "volume" in lines[5]
    assert "607.3" in lines[4]          # 636,834,220 B -> MiB
    assert "63.55" in lines[4]          # index build seconds
    assert "12,011,994" in lines[5]


def test_fleet3pt_table_keeps_not_run_explicit() -> None:
    rows = [_row(subtask="point_locate", copies=1, triangles=8_192_000,
                 seconds=None, throughput=None, status="not_run",
                 note="brute-force reference skipped")]
    body = render_fleet3pt_table(_fleet3pt_artifact(rows)).splitlines()[2]
    assert "not_run" in body
    assert "brute-force reference skipped" in body
    assert DASH in body


def test_fleet3pt_gates_parse_detail_and_mark_not_run() -> None:
    gates = [
        _gate("scaling_load", "exponent=1.17 (limit 1.50), R^2=0.999, 3 sizes 410,878..4,003,998",
              ratio=1.173),
        _gate("scaling_fleet_assemble[x2]",
              "not_run: 1 distinct size(s), need 3 to fit an exponent",
              passed=True, status="not_run"),
    ]
    table = render_fleet3pt_gates(_fleet3pt_artifact(gates=gates))
    lines = table.splitlines()
    load_row = [ln for ln in lines if "scaling_load" in ln][0]
    assert "1.17" in load_row and "1.50" in load_row
    assert "0.999" in load_row and "| 3 |" in load_row
    not_run_row = [ln for ln in lines if "fleet_assemble" in ln][0]
    assert "not_run" in not_run_row
    assert "pass" not in not_run_row.split("|")[2]


def test_committed_fleet3pt_artifact_measures_the_10m_band() -> None:
    """The anchor's missing piece: 10^6-10^7 band with *evaluable* gates."""
    artifact = load_artifact(FLEET3PT)
    biggest = max(r["triangles"] for r in artifact["results"])
    assert biggest >= 10_000_000, f"largest fleet rung is only {biggest:,}"
    scaling = [g for g in artifact["gates"] if g["name"].startswith("scaling_")]
    assert scaling, "no complexity gates in the fleet ladder"
    assert all(g["status"] == "measured" for g in scaling)
    assert artifact["gates_not_run"] == 0
    assert artifact["verdict"] == "pass"
    assert all(c["ok"] for c in artifact.get("volume_checks", []))


def test_committed_mesh_artifact_reaches_the_generator_limit() -> None:
    mesh = load_artifact(MESH)
    biggest = max(r["triangles"] for r in mesh["results"])
    assert biggest >= 8_000_000, f"largest single-body rung is only {biggest:,}"


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(pytest.main([__file__, "-q"]))