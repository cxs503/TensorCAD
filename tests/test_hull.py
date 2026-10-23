import json

import pytest

from tensorcad.hull import HullParameters, export_hull, generate_hull, validate_hull


def test_default_hull_is_deterministic_and_valid():
    first = generate_hull()
    second = generate_hull()
    assert first == second
    assert validate_hull(first) == []
    assert len(first.vertices) == 2 * 41 * 13
    assert len(first.faces) == 4 * 40 * 12
    assert len(first.offsets) == 41


def test_dimensions_change_geometry():
    baseline = generate_hull(HullParameters(length=10, beam=2, draft=0.8, nx=11, nz=5))
    wider = generate_hull(HullParameters(length=10, beam=3, draft=0.8, nx=11, nz=5))
    assert max(y for _, y, _ in wider.vertices) > max(y for _, y, _ in baseline.vertices)


@pytest.mark.parametrize("kwargs", [
    {"length": 0}, {"beam": -1}, {"draft": float("nan")}, {"nx": 2}, {"nz": 1}
])
def test_invalid_parameters_are_rejected(kwargs):
    values = {"length": 10, "beam": 2, "draft": 0.8, "nx": 11, "nz": 5}
    values.update(kwargs)
    with pytest.raises(ValueError):
        generate_hull(HullParameters(**values))


def test_export_creates_reproducible_artifacts(tmp_path):
    model = generate_hull(HullParameters(nx=9, nz=5))
    artifacts = export_hull(model, tmp_path)
    assert set(artifacts) == {"mesh", "offsets", "metadata"}
    assert all((tmp_path / name).is_file() for name in ("hull.obj", "offsets.csv", "metadata.json"))
    metadata = json.loads((tmp_path / "metadata.json").read_text())
    assert metadata["validation"] == []
    assert metadata["units"] == "m"
    assert metadata["vertex_count"] == 2 * 9 * 5
    assert (tmp_path / "hull.obj").read_text().count("\nv ") == 2 * 9 * 5
