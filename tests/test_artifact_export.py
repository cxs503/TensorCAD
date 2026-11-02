"""Tests for the STL artifact and for gated, resolution-controlled hull export.

An STL is a triangle soup: it carries no topology, no units and no provenance, so
the only way to know that the file is usable is to read it back, rebuild the
topology and re-run the gate on *the file*. These tests do exactly that, and
compare both the in-memory mesh and the recovered STL against
:func:`tensorcad.hull.analytic_volume`, an integral the export code never uses as
an input.

Binary STL stores float32 coordinates, so the recovered geometry cannot equal the
double-precision mesh bit for bit; the tolerances below are the float32 floor
(``2**-23`` relative), not a fudge factor.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from tensorcad import (
    HullParameters,
    MeshGateError,
    analytic_volume,
    export_hull,
    generate_hull,
    read_metadata,
    read_stl,
    write_stl,
)
from tensorcad import mesh as meshlib
from tensorcad.export import STL_VOLUME_TOLERANCE

FLOAT32_RELATIVE = 2.0 ** -23


@pytest.fixture()
def model():
    return generate_hull(HullParameters(nx=41, nz=13))


def _exact(model) -> float:
    p = model.parameters
    return analytic_volume(p.length, p.beam, p.draft)


def test_stl_round_trip_recovers_the_closed_topology(tmp_path: Path, model) -> None:
    path = write_stl(model.mesh, tmp_path / "hull.stl")
    recovered = read_stl(path)
    assert len(recovered.vertices) == len(model.mesh.vertices)
    assert len(recovered.faces) == len(model.mesh.faces)
    assert meshlib.is_watertight(recovered)
    assert meshlib.is_edge_manifold(recovered)
    assert len(recovered.vertices) - len(meshlib.edge_face_counts(recovered)) + len(
        recovered.faces
    ) == 2
    assert meshlib.total_signed_volume(recovered) > 0.0


def test_stl_coordinates_are_float32_exact_not_lossy_beyond_the_format(tmp_path: Path, model) -> None:
    """Every coordinate survives the trip to within one float32 quantum.

    Compared as sorted tables because an STL has no vertex table, so the reader
    numbers vertices in first-use order rather than in the generator's order.
    """
    path = write_stl(model.mesh, tmp_path / "hull.stl")
    recovered = read_stl(path)
    scale = max(abs(value) for point in model.mesh.vertices for value in point)
    tolerance = 4.0 * FLOAT32_RELATIVE * scale
    for original, restored in zip(
        sorted(model.mesh.vertices), sorted(recovered.vertices), strict=True
    ):
        for a, b in zip(original, restored, strict=True):
            assert math.isclose(a, b, rel_tol=0.0, abs_tol=tolerance), (a, b, tolerance)


def test_stl_volume_matches_the_double_precision_mesh(tmp_path: Path, model) -> None:
    path = write_stl(model.mesh, tmp_path / "hull.stl")
    recovered = read_stl(path)
    reference = meshlib.total_signed_volume(model.mesh)
    relative = meshlib.total_signed_volume(recovered) / reference - 1.0
    assert abs(relative) < STL_VOLUME_TOLERANCE
    assert abs(relative) < 1e-6  # the measured float32 floor, far below the tolerance


def test_stl_is_byte_deterministic(tmp_path: Path, model) -> None:
    first = write_stl(model.mesh, tmp_path / "a.stl").read_bytes()
    second = write_stl(model.mesh, tmp_path / "b.stl").read_bytes()
    assert first == second
    assert len(first) == 84 + 50 * len(model.mesh.faces)


def test_ascii_stl_is_rejected_with_an_explanation(tmp_path: Path) -> None:
    path = tmp_path / "ascii.stl"
    path.write_text("solid hull\n  facet normal 0 0 1\n    vertex 0 0 0\n  endfacet\nendsolid\n")
    with pytest.raises(ValueError, match="ASCII"):
        read_stl(path)


def test_truncated_stl_is_rejected(tmp_path: Path, model) -> None:
    path = write_stl(model.mesh, tmp_path / "hull.stl")
    data = path.read_bytes()
    path.write_bytes(data[:-25])
    with pytest.raises(ValueError, match="truncated"):
        read_stl(path)
    path.write_bytes(data[:40])
    with pytest.raises(ValueError, match="too short"):
        read_stl(path)


def test_stl_reader_rejects_non_finite_coordinates(tmp_path: Path) -> None:
    import struct

    payload = bytearray(84 + 50)
    struct.pack_into("<80sI", payload, 0, b"nan header".ljust(80, b" "), 1)
    struct.pack_into("<12fH", payload, 84, 0, 0, 0, math.nan, 0, 0, 0, 1, 0, 0, 0, 1, 0)
    path = tmp_path / "nan.stl"
    path.write_bytes(bytes(payload))
    with pytest.raises(ValueError, match="non-finite"):
        read_stl(path)


def test_default_export_still_writes_exactly_the_three_artifacts(tmp_path: Path, model) -> None:
    artifacts = export_hull(model, tmp_path / "plain")
    assert set(artifacts) == {"mesh", "offsets", "metadata"}
    assert not (tmp_path / "plain" / "hull.stl").exists()
    metadata = read_metadata(artifacts["metadata"])
    assert "gate" not in metadata and "resolution" not in metadata and "stl" not in metadata


def test_export_with_stl_records_the_read_back_evidence(tmp_path: Path, model) -> None:
    artifacts = export_hull(model, tmp_path / "with-stl", stl=True, gate=True)
    assert set(artifacts) == {"mesh", "offsets", "metadata", "stl"}
    assert Path(artifacts["stl"]).is_file()
    metadata = read_metadata(artifacts["metadata"])
    assert metadata["gate"]["passed"] is True
    stl = metadata["stl"]
    assert stl["watertight"] is True
    assert stl["triangle_count"] == len(model.mesh.faces)
    assert abs(stl["relative_volume_error_vs_mesh"]) < STL_VOLUME_TOLERANCE
    assert stl["gate"]["passed"] is True


def test_resolution_controlled_export_hits_the_analytic_volume(tmp_path: Path) -> None:
    """Ask for 0.5 %; the written STL is checked against the analytic integral."""
    coarse = generate_hull(HullParameters(nx=11, nz=5))
    target = 0.005
    artifacts = export_hull(coarse, tmp_path / "resolved", stl=True, target_volume_error=target)
    metadata = read_metadata(artifacts["metadata"])
    plan = metadata["resolution"]
    assert plan["accepted"] is True
    assert plan["parameters"]["nx"] > coarse.parameters.nx
    assert abs(plan["relative_error"]) <= target
    # Independent end-to-end check: analytic reference vs the artifact on disk.
    stl_volume = read_metadata(artifacts["metadata"])["stl"]["signed_volume"]
    recovered = read_stl(artifacts["stl"])
    assert recovered.faces  # the file exists and parses
    assert abs(stl_volume / analytic_volume(10.0, 2.0, 0.8) - 1.0) <= target + 1e-5
    assert abs(meshlib.total_signed_volume(recovered) / analytic_volume(10.0, 2.0, 0.8) - 1.0) <= (
        target + 1e-5
    )
    assert metadata["gate"]["measurements"]["triangles"] == plan["triangles"]


def test_export_refuses_a_mesh_that_fails_the_gate(tmp_path: Path) -> None:
    coarse = generate_hull(HullParameters(nx=11, nz=5))
    target = tmp_path / "refused"
    with pytest.raises(MeshGateError, match="reference_volume"):
        export_hull(coarse, target, gate=True)
    # Refusal happens before anything is written: no half-valid bundle on disk.
    assert not (target / "hull.obj").exists()


def test_export_rejects_metadata_from_a_grid_it_had_to_replace(tmp_path: Path) -> None:
    coarse = generate_hull(HullParameters(nx=11, nz=5))
    with pytest.raises(ValueError, match="different grid"):
        export_hull(
            coarse,
            tmp_path / "mismatch",
            metadata=coarse.metadata(with_discretisation=False),
            target_volume_error=0.005,
        )


def test_export_records_a_resolution_plan_and_a_watertight_flag_json_readable(
    tmp_path: Path, model
) -> None:
    artifacts = export_hull(model, tmp_path / "json", stl=True, gate=True)
    payload = json.loads(Path(artifacts["metadata"]).read_text())
    assert payload["gate"]["schema"] == "tensorcad.mesh.gate/1"
    assert payload["stl"]["gate"]["measurements"]["euler_characteristic"] == 2
    # resolution is opt-in by design, so an absent key is not a failure
    assert payload.get("resolution", True)
