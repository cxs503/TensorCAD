"""Subdivision-source rungs: import a *single connected* 10^7-triangle body.

``HullParameters`` caps one generated body at 8.19e6 triangles, so the 10^7 band
of GAP-CAD-D4-PERF-SCALE cannot be reached by the parametric generator alone and
the committed import ladder stopped at 4.0e6.  Uniform 1-to-4 midpoint
subdivision refines the *same discrete solid* (each planar parent triangle is
tiled exactly by its four children, midpoints welded by index identity), so it
gives one connected body at 13.2e6 triangles without changing shape.

These tests pin the four claims that make a refined rung a legitimate
measurement rather than a bigger number:

1. **Refinement does not change the solid.** Volume and topology are preserved
   per level, so a refined rung's volume error is the *parent's* discretisation
   error, not a new one.
2. **The round trip stays lossless at a refined rung.** The claim "importing
   10^7 triangles takes N seconds" is only valid if the parsed mesh is the same
   geometry that was written -- re-welding or deduplicating midpoints during the
   parse would silently change the triangle count.
3. **Refined rungs never share a baseline/glob key with the unrefined one.**
   A 4x larger mesh compared against the unrefined wall clock would invent a
   RED (or hide one), and a mixed-up OBJ file would report a parse time for the
   wrong triangle count.
4. **A refused rung is ``not_run`` with a reason.** Too many triangles for
   ``max_triangles``, or a sub-task that is not a function of the refined mesh,
   must say so instead of quietly measuring the unrefined mesh.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from tensorcad import perf_benchmark as pb
from tensorcad.hull import HullParameters, analytic_volume, generate_hull
from tensorcad.mesh import total_signed_volume
from tensorcad.mesh_subdivision import subdivide_mesh


@pytest.fixture()
def import_dir(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Redirect the OBJ scratch dir into tmp_path and drop every mesh cache."""
    target = tmp_path / "perf-import"
    target.mkdir()
    monkeypatch.setattr(pb, "IMPORT_DIR", target)
    pb.clear_mesh_cache()
    pb.clear_fleet_cache()
    pb.clear_subdiv_cache()
    yield target
    pb.clear_mesh_cache()
    pb.clear_fleet_cache()
    pb.clear_subdiv_cache()


def _hull(nx: int = 41, nz: int = 13):
    return generate_hull(
        HullParameters(length=10.0, beam=2.0, draft=0.8, nx=nx, nz=nz)
    ).mesh


def test_refinement_preserves_the_solid_volume_and_topology() -> None:
    mesh = _hull()
    volume = total_signed_volume(mesh)
    reference = analytic_volume(10.0, 2.0, 0.8)

    refined = subdivide_mesh(mesh, 2)
    assert len(refined.faces) == len(mesh.faces) * 4 ** 2
    # Exact tiling of planar parents: the discrete solid may not change.
    assert total_signed_volume(refined) == pytest.approx(volume, rel=1e-14)
    # ... and the inherited sampling error is the parent's, unchanged.
    assert (volume - reference) / reference == pytest.approx(
        (total_signed_volume(refined) - reference) / reference, rel=1e-12
    )


def test_subdivided_accessor_is_cached_and_identity_at_level_zero(import_dir: Path) -> None:
    base = pb._subdivided_for(41, 13, 1, 0)
    assert base is pb._mesh_for(41, 13, 1)
    level2 = pb._subdivided_for(41, 13, 1, 2)
    assert level2 is pb._subdivided_for(41, 13, 1, 2)  # cached, not rebuilt
    assert len(level2.faces) == len(base.faces) * 16
    # Incremental build: level 2 reuses level 1 rather than re-refining the base.
    pb._subdivided_for(41, 13, 1, 1)
    level3 = pb._subdivided_for(41, 13, 1, 3)
    assert len(level3.faces) == len(base.faces) * 64


def test_obj_round_trip_is_lossless_on_a_refined_rung(import_dir: Path) -> None:
    write = pb._run_obj_write(21, 9, levels=2)
    parse = pb._run_obj_read(21, 9, levels=2)

    assert write.status == pb.MEASURED and parse.status == pb.MEASURED
    assert write.levels == 2 and parse.levels == 2
    assert write.triangles == parse.triangles == 678 * 16
    assert pb.import_obj_path(21, 9, levels=2).name.endswith("_l2.obj")
    # A parse that re-welded midpoints would change the count; a parse that lost
    # precision would move the volume. Both are checked by the sub-task itself.
    assert parse.volume_rel_err == pytest.approx(write_level0_error(import_dir), rel=1e-6)


def write_level0_error(import_dir: Path) -> float:
    """Volume error of the *unrefined* mesh at the same base resolution."""
    result = pb._run_obj_read(21, 9, levels=0)
    return result.volume_rel_err


def test_refined_rung_keeps_its_own_baseline_key(import_dir: Path, tmp_path: Path) -> None:
    assert pb.measurement_key("obj_read", 41, 13) == "obj_read:41x13"
    assert pb.measurement_key("obj_read", 41, 13, levels=2) == "obj_read:41x13:l2"

    result = pb.run_perf_benchmark(
        ladder=((21, 9),), voxel_ladder=(), subtasks=("obj_read",),
        n_queries=8, subdivide_levels=(0, 1),
    )
    baseline = tmp_path / "baseline.json"
    baseline.write_text(json.dumps({
        "schema": pb.BENCHMARK_SCHEMA,
        "environment": result.environment,
        "results": [
            {"subtask": "obj_read", "nx": 21, "nz": 9, "copies": 1,
             "seconds": 1e-6}  # absurdly small: any comparison would be RED
        ],
    }), encoding="utf-8")

    gated, failed = pb.evaluate_gates(result, baseline, regression_factor=2.0)
    names = [g.name for g in gated.gates]
    # The unrefined rung is compared (and is RED against the 1 µs baseline);
    # the refined rung has no baseline and must not inherit the unrefined time.
    assert "regression_obj_read:21x9" in names
    assert "regression_obj_read:21x9:l1" not in names


def test_refined_rung_beyond_max_triangles_is_not_run(import_dir: Path) -> None:
    result = pb.run_perf_benchmark(
        ladder=((41, 13),), voxel_ladder=(), subtasks=("obj_write",),
        n_queries=8, subdivide_levels=(0, 3), max_triangles=100_000,
    )
    by_level = {r.levels: r for r in result.results}
    assert by_level[0].status == pb.MEASURED
    refused = by_level[3]
    assert refused.status == pb.NOT_RUN
    assert refused.seconds is None
    assert "max_triangles" in refused.note
    assert refused.triangles == 4 * 40 * 12 * 64  # the size it refused to run


def test_task_that_ignores_the_refined_mesh_is_not_run(import_dir: Path) -> None:
    result = pb.run_perf_benchmark(
        ladder=((41, 13),), voxel_ladder=(), subtasks=("load", "volume"),
        n_queries=8, subdivide_levels=(1,),
    )
    by_task = {r.subtask: r for r in result.results}
    # `load` regenerates its own source, so it cannot be labelled "refined".
    assert by_task["load"].status == pb.NOT_RUN
    assert "subdivide_levels applies only to" in by_task["load"].note
    assert by_task["volume"].status == pb.MEASURED
    assert by_task["volume"].levels == 1


def test_scaling_gates_are_measured_across_levels(import_dir: Path) -> None:
    result = pb.run_perf_benchmark(
        ladder=((21, 9),), voxel_ladder=(), subtasks=("obj_write", "obj_read", "volume"),
        n_queries=8, subdivide_levels=(0, 1, 2),
    )
    gates = {g.name: g for g in result.gates}
    for subtask in ("obj_write", "obj_read", "volume"):
        gate = gates[f"scaling_{subtask}"]
        assert gate.status == pb.MEASURED, gate.detail
        # 678 -> 10,848 triangles: three distinct sizes, O(N) path.
        assert float(gate.ratio) < pb.SCALING_LIMITS[subtask], gate.detail
        assert "3 sizes" in gate.detail


def test_clear_import_files_removes_refined_files_too(import_dir: Path) -> None:
    from tensorcad.export import write_obj

    # Written directly (not via _run_obj_write, which prunes the other rung):
    # both rungs must be visible to the cleanup glob at the same time.
    write_obj(pb._subdivided_for(21, 9, 1, 0), pb.import_obj_path(21, 9, levels=0))
    write_obj(pb._subdivided_for(21, 9, 1, 2), pb.import_obj_path(21, 9, levels=2))
    removed = pb.clear_import_files()
    assert sorted(removed) == ["hull_21x9_b1.obj", "hull_21x9_b1_l2.obj"]
    assert list(import_dir.glob("*.obj")) == []