"""OBJ file import/export sub-tasks of the large-geometry benchmark.

These tests pin three properties that a "file import is fast" claim would
otherwise hide:

1. **The round trip is lossless.** ``write_obj`` writes 17 significant digits,
   so parsing must return bit-identical coordinates and the same triangle count;
   the volume therefore equals the in-memory mesh's volume exactly.
2. **The benchmark measures what it says.** ``obj_write``/``obj_read`` come back
   ``measured`` with a throughput and a counted file size, and ``obj_read``
   carries the as-imported volume error against the analytic reference so a
   fast-but-wrong parse cannot be reported as a good number.
3. **Failure is loud.** A file that lost triangles raises instead of returning a
   small ``seconds`` value; a filesystem with no room returns ``not_run`` with a
   reason (never a fake pass), and the 2x regression gate is RED when the parse
   really did double in time.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

from tensorcad import perf_benchmark as pb
from tensorcad.export import read_obj, write_obj
from tensorcad.hull import HullParameters, analytic_volume, generate_hull
from tensorcad.mesh import total_signed_volume


@pytest.fixture()
def import_dir(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Redirect the benchmark's OBJ scratch directory into the test's tmp_path."""
    target = tmp_path / "perf-import"
    target.mkdir()
    monkeypatch.setattr(pb, "IMPORT_DIR", target)
    pb.clear_mesh_cache()
    yield target
    pb.clear_mesh_cache()


def _mesh(nx: int = 41, nz: int = 13):
    return generate_hull(HullParameters(length=10.0, beam=2.0, draft=0.8, nx=nx, nz=nz)).mesh


def test_obj_round_trip_is_lossless(import_dir: Path) -> None:
    mesh = _mesh()
    path = write_obj(mesh, import_dir / "round_trip.obj")
    parsed = read_obj(path)

    assert len(parsed.faces) == len(mesh.faces)
    assert len(parsed.vertices) == len(mesh.vertices)
    # 17 significant digits: the coordinate round trip is exact, so the volume
    # is bit-identical rather than "within tolerance".
    assert parsed.vertices == mesh.vertices
    assert parsed.faces == mesh.faces
    assert total_signed_volume(parsed) == total_signed_volume(mesh)


def test_obj_subtasks_are_measured_with_throughput(import_dir: Path) -> None:
    written = pb._run_obj_write(41, 13)
    read = pb._run_obj_read(41, 13)

    assert written.status == pb.MEASURED and read.status == pb.MEASURED
    assert written.seconds > 0.0 and read.seconds > 0.0
    assert written.throughput and written.throughput > 0.0
    assert read.throughput and read.throughput > 0.0
    assert written.triangles == read.triangles == len(_mesh().faces)
    assert "MiB" in written.note and "B/tri" in written.note
    # The as-imported geometry is judged against the analytic reference, not
    # against the file: that is what makes the row usable in an error table.
    assert "rel_err" in read.note and "analytic" in read.note
    assert pb.import_obj_path(41, 13).exists()


def test_obj_read_prepares_the_file_when_missing_and_says_so(import_dir: Path) -> None:
    assert not pb.import_obj_path(81, 21).exists()
    result = pb._run_obj_read(81, 21)
    assert result.status == pb.MEASURED
    assert "prepared in" in result.note  # preparation time is disclosed, untimed


def test_obj_read_raises_when_the_file_lost_triangles(import_dir: Path) -> None:
    mesh = _mesh()
    path = pb.import_obj_path(41, 13)
    write_obj(mesh, path)
    # Keep the syntax valid but drop the last face: exactly the silent-loss case
    # a bare wall-clock number would hide.
    lines = path.read_text(encoding="utf-8").splitlines()
    face_indices = [i for i, line in enumerate(lines) if line.startswith("f ")]
    path.write_text("\n".join(lines[: face_indices[-1]]) + "\n", encoding="utf-8")

    with pytest.raises(ValueError, match="round trip lost geometry"):
        pb._run_obj_read(41, 13)


def test_obj_write_is_not_run_when_the_disk_is_full(import_dir: Path, monkeypatch) -> None:
    monkeypatch.setattr(pb, "_free_mb", lambda path: 0.0)
    result = pb._run_obj_write(41, 13)
    assert result.status == pb.NOT_RUN
    assert result.seconds is None
    assert "not enough free space" in result.note
    assert not pb.import_obj_path(41, 13).exists()


def test_import_subtasks_get_measured_scaling_gates(import_dir: Path) -> None:
    ladder = ((41, 13), (81, 21), (161, 41))
    result = pb.run_perf_benchmark(
        ladder=ladder,
        voxel_ladder=(),
        subtasks=("load", "obj_write", "obj_read"),
        n_queries=16,
    )
    names = {gate.name: gate for gate in result.gates}
    for subtask in ("obj_write", "obj_read"):
        gate = names[f"scaling_{subtask}"]
        assert gate.status == pb.MEASURED, gate.detail
        # O(N) serialisation/parse: an exponent near 1, far from O(N^2).
        exponent = float(re.search(r"exponent=([-\d.]+)", gate.detail).group(1))
        assert 0.0 < exponent < pb.SCALING_LIMITS[subtask], gate.detail


def test_regression_gate_flags_a_doubled_parse(import_dir: Path, tmp_path: Path) -> None:
    result = pb.run_perf_benchmark(
        ladder=((41, 13),), voxel_ladder=(), subtasks=("obj_write", "obj_read"),
        n_queries=16,
    )
    baseline = tmp_path / "baseline.json"
    payload = {
        "schema": pb.BENCHMARK_SCHEMA,
        "environment": result.environment,
        "results": [
            {"subtask": r.subtask, "nx": r.nx, "nz": r.nz, "copies": r.copies,
             "seconds": r.seconds / 2.0}  # baseline is half today's time -> 2x
            for r in result.results
        ],
    }
    baseline.write_text(json.dumps(payload), encoding="utf-8")

    gated, failed = pb.evaluate_gates(result, baseline, regression_factor=2.0)
    obj_read_gate = next(g for g in gated.gates if g.name == "regression_obj_read:41x13")
    assert obj_read_gate.ratio == pytest.approx(2.0, rel=1e-6)
    assert obj_read_gate.passed is True  # ratio == limit is still within budget
    assert failed == []

    stricter, failed_strict = pb.evaluate_gates(result, baseline, regression_factor=1.5)
    names = [g.name for g in strict_only(stricter)]
    assert "regression_obj_read:41x13" in names
    assert "regression_obj_read:41x13" in failed_strict


def strict_only(result) -> list:
    by_name = {gate.name: gate for gate in result.gates}
    return [by_name["regression_obj_read:41x13"]]


def test_clear_import_files_removes_only_benchmark_files(import_dir: Path) -> None:
    write_obj(_mesh(), pb.import_obj_path(41, 13))
    keep = import_dir / "handwritten.obj"
    keep.write_text("# not ours\n", encoding="utf-8")
    removed = pb.clear_import_files()
    assert removed == ["hull_41x13_b1.obj"]
    assert keep.exists()
    assert not pb.import_obj_path(41, 13).exists()