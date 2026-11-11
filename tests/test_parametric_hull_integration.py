"""Integration tests: parametric hull -> FEM/CFD handoff, and the benchmark record.

Three things are pinned here that the unit tests cannot see:

1. **Orientation.** The parametric grid used to wind its triangles inward (the closed
   Wigley body measured -2.7734 m^3 against +2.7734 outward), which inverted every load
   sign in a consumer and made the FEM adapter reject the surface.  The sign is now
   fixed at the source, so both the closed body and the hull patch alone must integrate
   positive.
2. **The handoff.** :func:`tensorcad.adapters.fem.parametric_hull_surface` is the call a
   solver-side consumer makes; it must produce a watertight, outward, sphere-topology
   surface whose volume equals the hydrostatics of the same mesh.
3. **The published benchmark record.** ``benchmarks/verified/parametric_hull/result.json``
   is a committed artifact, so a test has to fail when it drifts from the record it was
   generated from -- otherwise the numbers in the repository can silently rot.
"""
from __future__ import annotations

import json
import pathlib

import pytest

from tensorcad import mesh as meshlib
from tensorcad.adapters.fem import parametric_hull_surface, read_gmsh, write_gmsh
from tensorcad.parametric_hull import (
    build_parametric_hull,
    wigley_analytic_properties,
    wigley_cross_sections,
)

ROOT = pathlib.Path(__file__).resolve().parents[1]
L, B, T = 10.0, 1.0, 0.625


@pytest.fixture(scope="module")
def hull():
    sections = wigley_cross_sections(L, B, T, 11, 12)
    return build_parametric_hull(
        sections, mesh_n_u=97, mesh_n_v=16, close_waterplane=True
    )


def test_closed_hull_surface_is_oriented_outward(hull):
    volume = meshlib.total_signed_volume(hull.mesh_data.mesh)
    assert volume > 0.0
    # Same body, measured two ways: divergence volume on the mesh, hydrostatics on the
    # clipped wetted surface.  They agree because the lid sits on the waterline.
    assert volume == pytest.approx(hull.hydrostatics.displaced_volume, rel=1e-12)


def test_hull_patch_alone_is_oriented_outward(hull):
    volume = meshlib.total_signed_volume(hull.mesh_data.hull_mesh)
    assert volume > 0.0


def test_parametric_hull_surface_is_watertight_and_carries_quality(hull):
    surface = parametric_hull_surface(hull)
    assert surface.watertight is True
    assert surface.edge_manifold is True
    assert surface.euler_characteristic == 2
    analytic = wigley_analytic_properties(L, B, T)
    assert surface.volume_m3 == pytest.approx(analytic["displaced_volume"], rel=5e-3)
    assert surface.volume_m3 > 0.0  # outward: the adapter would have raised otherwise
    recorded = surface.quality["parametric_hull"]
    assert recorded["n_u"] == hull.mesh_data.n_u
    assert recorded["n_v"] == hull.mesh_data.n_v
    assert recorded["mesh_quality"]["min_scaled_jacobian"] > 0.0
    assert "waterline" in surface.coordinate_frame


def test_parametric_hull_surface_rejects_open_shell():
    sections = wigley_cross_sections(L, B, T, 11, 12)
    open_hull = build_parametric_hull(
        sections, mesh_n_u=41, mesh_n_v=10, close_waterplane=False
    )
    with pytest.raises(ValueError, match="close_waterplane"):
        parametric_hull_surface(open_hull)


def test_parametric_hull_surface_rejects_wrong_type():
    with pytest.raises(TypeError, match="ParametricHull"):
        parametric_hull_surface(meshlib.weld)  # any non-hull object


def test_gmsh_roundtrip_preserves_the_mesh(hull, tmp_path):
    surface = parametric_hull_surface(hull)
    path = write_gmsh(surface, tmp_path / "hull.msh")
    parsed = read_gmsh(path)
    assert parsed["nodes"] == len(surface.nodes)
    assert parsed["triangles"] == len(surface.triangles)
    reread = read_gmsh(path)
    assert reread["physical_names"]  # one tag per region survived the write


# --------------------------------------------------------------------------- #
# committed benchmark record
# --------------------------------------------------------------------------- #
BENCH = ROOT / "benchmarks" / "verified" / "parametric_hull"


def _load(name: str) -> dict:
    path = BENCH / name
    if not path.is_file():
        pytest.skip(f"{path.relative_to(ROOT)} not generated yet; run the benchmark")
    return json.loads(path.read_text())


def test_benchmark_record_agrees_with_verification_record():
    result = _load("result.json")
    verification = _load("verification.json")
    by_id = {case["id"]: case for case in result["cases"]}
    assert by_id["b_wigley_displacement"]["measured"] == pytest.approx(
        verification["criterion_b"]["measured_volume"], rel=1e-12
    )
    assert by_id["c_wigley_waterplane"]["measured"] == pytest.approx(
        verification["criterion_c"]["measured_waterplane_area"], rel=1e-12
    )
    assert by_id["d_suboff_offsets"]["measured"] == pytest.approx(
        verification["criterion_d"]["value"], rel=1e-12
    )
    assert by_id["e_min_scaled_jacobian"]["measured"] == pytest.approx(
        verification["criterion_e"]["value"]["min_scaled_jacobian"], rel=1e-12
    )


def test_benchmark_record_is_a_pass_and_cites_its_repro_command():
    result = _load("result.json")
    assert result["all_passed"] is True
    assert result["repro_cmd"].endswith("run_benchmark.py")
    assert all(case["passed"] for case in result["cases"])
    # Criterion (b) and (c) thresholds from the node text, in the record.
    by_id = {case["id"]: case for case in result["cases"]}
    assert by_id["b_wigley_displacement"]["threshold"] == pytest.approx(5e-3)
    assert by_id["c_wigley_waterplane"]["threshold"] == pytest.approx(1e-2)
    assert abs(by_id["b_wigley_displacement"]["error"]) < 5e-3
    assert abs(by_id["c_wigley_waterplane"]["error"]) < 1e-2


def test_benchmark_figures_exist():
    result = _load("result.json")
    figures = result["artifacts"]["figures"]
    assert len(figures) >= 6
    for figure in figures:
        assert (ROOT / figure["path"]).is_file(), figure["path"]
