"""End-to-end tests for the neural geometry pipeline (requires the ``neural`` extra)."""
from __future__ import annotations

import json
import shutil
from pathlib import Path

import pytest

torch = pytest.importorskip("torch", reason="the neural geometry runtime needs the 'neural' extra")

from tensorcad.neural.objects import (  # noqa: E402
    NeuralObjectError,
    Tolerances,
    apply_dent_edit,
    extract_object,
    fit_sphere,
)
from tensorcad.neural_geometry import ObjectStore, load_bundle  # noqa: E402

LOOSE = Tolerances(
    field_mae=1e9,
    boundary_error_mean=1e9,
    ray_coverage=0.0,
    dent_peak_error=1e9,
    dent_region_error_mean=1e9,
    protected_mean_displacement=1e9,
    protected_p99_displacement=1e9,
    protected_max_displacement=1e9,
)

STEPS = 200
RAYS = 96


@pytest.fixture(scope="module")
def base_bundle(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Fit one object for the whole module; each test works on its own copy."""
    directory = tmp_path_factory.mktemp("neural-base")
    target = directory / "sphere.ngo"
    fit_sphere(target, radius=0.5, steps=250, seed=7, tolerance=LOOSE, rays=RAYS)
    return target


@pytest.fixture()
def fitted(tmp_path: Path, base_bundle: Path) -> Path:
    target = tmp_path / "sphere.ngo"
    shutil.copytree(base_bundle, target)
    return target


def test_fit_publishes_revision_one_with_metadata(fitted: Path) -> None:
    target = fitted
    manifest = load_bundle(target)
    assert manifest.revision == 1
    assert manifest.lineage.parent_revision is None
    assert manifest.validation_status == "passed"
    assert manifest.units == "m"
    assert manifest.field["quantity"] == "truncated_signed_distance"
    store = ObjectStore(target)
    assert [entry["status"] for entry in store.history()] == ["committed"]
    assert store.active_revision() == 1


def test_fit_refuses_to_overwrite_an_existing_object(fitted: Path) -> None:
    with pytest.raises(NeuralObjectError, match="already holds revision"):
        fit_sphere(fitted, steps=10, tolerance=LOOSE, rays=8)


def test_edit_commits_a_new_revision_and_records_lineage(fitted: Path) -> None:
    target = fitted
    parent_hash = ObjectStore(target).active_weights_sha256()
    result = apply_dent_edit(
        target, depth=0.05, sigma=0.25, steps=STEPS, seed=9, tolerance=LOOSE, rays=RAYS
    )
    assert result.committed
    assert result.parent_revision == 1
    assert result.active_revision == 2
    manifest = load_bundle(target)
    assert manifest.revision == 2
    assert manifest.lineage.parent_revision == 1
    assert manifest.lineage.root_revision == 1
    assert manifest.raw["training"]["parent_weights_sha256"] == parent_hash
    entries = ObjectStore(target).history()
    assert [entry["revision"] for entry in entries] == [1, 2]


def test_failed_edit_rolls_back_and_leaves_the_parent_untouched(fitted: Path) -> None:
    target = fitted
    store = ObjectStore(target)
    before = store.active_weights_sha256()
    before_manifest = (store.root / "manifest.json").read_bytes()

    result = apply_dent_edit(
        target, depth=0.30, sigma=0.08, steps=STEPS, seed=11, tolerance=Tolerances(), rays=RAYS
    )
    assert not result.committed
    assert result.status == "rolled_back"
    assert result.active_revision == 1
    assert result.failed_checks, "a rejection must name the gates it broke"
    assert store.active_weights_sha256() == before
    assert (store.root / "manifest.json").read_bytes() == before_manifest
    assert store.active_revision() == 1

    entries = store.history()
    assert [entry["status"] for entry in entries] == ["committed", "rolled_back"]
    assert entries[-1]["failed_checks"] == list(result.failed_checks)
    assert entries[-1]["revision"] == 2, "the attempt is recorded against the revision it tried to publish"


def test_rolled_back_metrics_stay_json_safe(fitted: Path) -> None:
    target = fitted
    apply_dent_edit(target, depth=0.45, sigma=0.06, steps=STEPS, seed=13, tolerance=Tolerances(), rays=RAYS)
    text = (target / "history.json").read_text(encoding="utf-8")
    json.loads(text, parse_constant=lambda token: pytest.fail(f"non-standard JSON token {token}"))
    entry = json.loads(text)[-1]
    assert entry["status"] == "rolled_back"
    assert isinstance(entry["metrics"], dict)


def test_edit_can_attach_derived_meshes_at_commit(fitted: Path) -> None:
    target = fitted
    result = apply_dent_edit(
        target,
        depth=0.05,
        sigma=0.25,
        steps=STEPS,
        seed=17,
        tolerance=LOOSE,
        rays=RAYS,
        extract_resolution=20,
    )
    assert result.committed
    manifest = load_bundle(target)
    kinds = {artifact["kind"] for artifact in manifest.derived}
    formats = {artifact["format"] for artifact in manifest.derived}
    assert kinds == {"surface_mesh"}
    assert formats == {"obj", "stl"}
    for artifact in manifest.derived:
        path = target / artifact["path"]
        assert path.is_file()
        assert path.stat().st_size > 0
        assert artifact["metrics"]["watertight"] is True


def test_extract_from_the_active_revision_writes_closed_meshes(fitted: Path) -> None:
    target = fitted
    mesh, report, records = extract_object(target, resolution=24, formats=("obj", "stl"))
    assert report.watertight
    assert report.edge_manifold
    assert report.euler_characteristic == 2
    assert report.signed_volume > 0.0
    assert len(records) == 2
    for record in records:
        assert (target / record["path"]).is_file()
        assert record["sha256"]
    assert mesh.units == "m"
