"""Tests for artifact writing, reading, and the reproducibility record."""
from __future__ import annotations

import itertools
import json
import math
from pathlib import Path

import pytest

from tensorcad import HullParameters, export_hull, generate_hull, read_metadata, read_obj, write_obj
from tensorcad import mesh as meshlib
from tensorcad.export import (
    ROUND_TRIP_PRECISION,
    ExportArtifacts,
    _format,
    write_metadata,
    write_offsets_csv,
)
from tensorcad.geometry import Mesh


@pytest.fixture()
def model():
    return generate_hull(HullParameters(nx=11, nz=5))


def test_obj_round_trips_to_the_identical_model(tmp_path: Path, model) -> None:
    """Coordinates must survive writing and reading bit-for-bit.

    The closed topology relies on vertices shared by *exact* position, so a file
    that only approximately reproduces them produces a different -- and open --
    model. 17 significant digits is what that requires, and this test is the
    reason it is the default rather than the usual 6.
    """
    path = write_obj(model.mesh, tmp_path / "hull.obj")
    reloaded = read_obj(path)
    assert reloaded.vertices == model.vertices
    assert reloaded.faces == model.faces
    assert meshlib.is_watertight(reloaded)


def test_reduced_precision_is_documented_as_lossy(tmp_path: Path, model) -> None:
    exact = read_obj(write_obj(model.mesh, tmp_path / "exact.obj", precision=ROUND_TRIP_PRECISION))
    lossy = read_obj(write_obj(model.mesh, tmp_path / "lossy.obj", precision=4))
    assert len(lossy.vertices) == len(exact.vertices)
    assert lossy.vertices != exact.vertices, "4 digits perturbs the model"
    assert lossy.faces == exact.faces, "index references survive; only coordinates drift"
    drift = max(
        abs(a - b)
        for pa, pb in zip(exact.vertices, lossy.vertices, strict=False)
        for a, b in zip(pa, pb, strict=False)
    )
    assert drift < 1e-3


def test_obj_header_states_units_and_frame(tmp_path: Path, model) -> None:
    text = (write_obj(model.mesh, tmp_path / "h.obj")).read_text()
    head = text.splitlines()[:6]
    assert any("# units: m" in line for line in head)
    assert any("coordinate frame" in line for line in head)
    assert "triangles:" in head[3]
    with pytest.raises(ValueError, match="precision"):
        write_obj(model.mesh, tmp_path / "x.obj", precision=0)


def test_negative_zero_never_reaches_a_file(tmp_path: Path) -> None:
    mesh = Mesh(vertices=((-0.0, -0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)), faces=((0, 1, 2),))
    line = next(l for l in write_obj(mesh, tmp_path / "z.obj").read_text().splitlines() if l.startswith("v "))
    assert line.split() == ["v", "0", "0", "0"], "no '-0' token for a reader to choke on"
    assert _format(-0.0, 17) == "0"


def test_read_obj_handles_quads_and_negative_indices(tmp_path: Path) -> None:
    source = tmp_path / "quad.obj"
    source.write_text(
        "# comment\nv 0 0 0\nv 1 0 0\nv 1 1 0\nv 0 1 0\nvn 0 0 1\nf 4 3 2 1\nf -1 -2 -3\n",
        encoding="utf-8",
    )
    mesh = read_obj(source)
    assert len(mesh.vertices) == 4
    assert mesh.faces == ((3, 2, 1), (3, 1, 0), (3, 2, 1))
    for bad in ("v 1 2\n", "f 1 2\n"):
        broken = tmp_path / "broken.obj"
        broken.write_text(f"v 0 0 0\nv 1 0 0\nv 1 1 0\n{bad}", encoding="utf-8")
        with pytest.raises(ValueError):
            read_obj(broken)


def test_offsets_csv_is_machine_readable(tmp_path: Path, model) -> None:
    path = write_offsets_csv(model.offsets, tmp_path / "offsets.csv")
    rows = path.read_text().strip().splitlines()
    assert rows[0] == "x_m,overall_breadth_m"
    assert len(rows) == 1 + model.parameters.nx
    parsed = [[float(value) for value in row.split(",")] for row in rows[1:]]
    assert parsed[0][1] == 0.0 and parsed[-1][1] == 0.0
    assert max(row[1] for row in parsed) == pytest.approx(model.parameters.beam)
    assert all(a[0] < b[0] for a, b in itertools.pairwise(parsed))


def test_metadata_never_emits_invalid_json(tmp_path: Path) -> None:
    """NaN is not JSON; it must surface as null rather than break other parsers."""
    payload = {"volume": float("nan"), "nested": {"inf": math.inf, "list": [1.0, float("nan")]}}
    path = write_metadata(payload, tmp_path / "m.json")
    loaded = json.loads(path.read_text())
    assert loaded["nested"]["list"][1] is None
    assert loaded["nested"]["inf"] is None
    text = path.read_text()
    assert "NaN" not in text and "Infinity" not in text
    assert text.endswith("\n")
    assert list(loaded) == sorted(loaded), "keys are sorted so two runs diff cleanly"


def test_export_writes_the_three_artifacts(tmp_path: Path, model) -> None:
    artifacts = export_hull(model, tmp_path / "run")
    assert ExportArtifacts(**artifacts).as_dict() == artifacts
    for name in ("hull.obj", "offsets.csv", "metadata.json"):
        assert (tmp_path / "run" / name).is_file()
    metadata = read_metadata(artifacts["metadata"])
    assert metadata["validation"]["passed"] is True
    assert metadata["mesh"]["vertex_count"] == len(model.vertices)
    assert metadata["schema_version"] == "2.0"
    assert "discretisation" not in metadata, "convergence reporting is opt-in"


def test_export_with_convergence_records_measured_error(tmp_path: Path, model) -> None:
    artifacts = export_hull(model, tmp_path / "with-study", convergence=True)
    metadata = read_metadata(artifacts["metadata"])
    study = metadata["discretisation"]
    assert study["status"] == "measured"
    assert study["refined_grid"] == {"nx": 21, "nz": 9}
    assert 0.0 < study["relative_change"]["displacement_volume"] < 0.2
    assert Path(artifacts["mesh"]).read_text() == (tmp_path / "with-study" / "hull.obj").read_text()
