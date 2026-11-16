"""Mesh-derived voxelisation ladder: artifact contract, gates, doc rendering."""
from __future__ import annotations

import json

import pytest

from tensorcad.perf_voxel_mesh import (
    SCHEMA,
    check_section,
    render_cell_table,
    render_section,
    render_triangle_table,
    run_voxelize_mesh_benchmark,
    scaling_gates,
)


@pytest.fixture(scope="module")
def artifact() -> dict:
    return run_voxelize_mesh_benchmark(
        bases=((21, 9),),
        levels_ladder=(0, 1),
        grids=((8, 4, 2),),
        band_factor=1.0,
    )


def test_artifact_schema_and_measured_rows(artifact: dict) -> None:
    assert artifact["schema"] == SCHEMA
    measured = [r for r in artifact["results"] if r["status"] == "measured"]
    assert len(measured) == 2
    sizes = sorted(r["triangles"] for r in measured)
    assert sizes[0] < sizes[-1]
    for row in measured:
        assert row["seconds"] > 0
        assert row["grid_cells"] == 8 * 4 * 2
        assert row["band"] > 0
        assert row["voxel_volume"] > 0
        assert row["volume_rel_err"] is not None
        assert row["us_per_query"] == pytest.approx(1e6 * row["seconds"] / row["grid_cells"], rel=0.01)


def test_only_two_sizes_means_not_run_not_pass(artifact: dict) -> None:
    names = {g["name"]: g for g in artifact["gates"]}
    assert names["scaling_voxelize_mesh"]["status"] == "not_run"
    assert names["scaling_voxelize_mesh"]["passed"] is None
    assert names["voxel_volume_accuracy"]["status"] == "measured"
    assert names["scaling_voxelize_mesh_cells"]["status"] == "not_run"


def test_quadratic_cost_goes_red() -> None:
    rows = [
        {"grid": [8, 4, 2], "triangles": n, "seconds": (n / 1000.0) ** 2,
         "status": "measured", "grid_cells": 64, "volume_rel_err": 0.001}
        for n in (1000, 2000, 4000)
    ]
    gates = {g["name"]: g for g in scaling_gates(rows)}
    assert gates["scaling_voxelize_mesh"]["status"] == "measured"
    assert gates["scaling_voxelize_mesh"]["passed"] is False


def test_linear_cost_passes_and_is_fitted() -> None:
    rows = [
        {"grid": [8, 4, 2], "triangles": n, "seconds": n / 1000.0,
         "status": "measured", "grid_cells": 64, "volume_rel_err": 0.001}
        for n in (1000, 2000, 4000)
    ]
    gates = {g["name"]: g for g in scaling_gates(rows)}
    assert gates["scaling_voxelize_mesh"]["passed"] is True
    assert gates["scaling_voxelize_mesh"]["ratio"] == pytest.approx(1.0, abs=1e-6)


def test_volume_error_above_3pct_goes_red() -> None:
    rows = [{
        "grid": [8, 4, 2], "triangles": 1000, "seconds": 1.0, "status": "measured",
        "grid_cells": 64, "volume_rel_err": 0.05,
    }]
    gates = {g["name"]: g for g in scaling_gates(rows)}
    assert gates["voxel_volume_accuracy"]["passed"] is False


def test_section_renders_and_check_detects_drift(artifact: dict) -> None:
    section = render_section(artifact)
    assert render_triangle_table(artifact) in section
    assert render_cell_table(artifact) in section
    doc = f"intro\n<!-- perf-meshvoxel:section begin -->\n{section}\n<!-- perf-meshvoxel:section end -->\ntail\n"
    assert check_section(doc, artifact) == []
    drifted = doc.replace("| 三角面 |", "| 三角面数 |")
    assert check_section(drifted, artifact)


def test_band_factor_must_be_positive() -> None:
    with pytest.raises(ValueError):
        run_voxelize_mesh_benchmark(bases=((21, 9),), levels_ladder=(0,), grids=((8, 4, 2),), band_factor=0.0)


def test_rung_above_max_triangles_is_not_run() -> None:
    artifact = run_voxelize_mesh_benchmark(
        bases=((21, 9),), levels_ladder=(0, 1), grids=((8, 4, 2),), max_triangles=2000,
    )
    statuses = [r["status"] for r in artifact["results"]]
    assert "not_run" in statuses
    not_run = [r for r in artifact["results"] if r["status"] == "not_run"][0]
    assert "max_triangles" in not_run["note"]
    assert json.dumps(artifact)  # JSON-serialisable evidence record
