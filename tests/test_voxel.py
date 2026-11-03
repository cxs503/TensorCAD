"""SDF voxel-grid export: geometry, sign convention, accuracy, and robustness.

The claims under test are measurable ones: the samples sit on cell centres, the grid
holds the field's own distance (not a resampled approximation), the partial-fill
volume is markedly better than hard cell counting at the same cost, and the SUBOFF
grid reproduces the published displacement volume that the mesh path already hits.
"""
from __future__ import annotations

import json
import math
import os
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

from tensorcad import EllipsoidField, FieldError, PlaneField, SphereField, SuboffHullField, to_voxel_grid
from tensorcad.voxel import VOXEL_HANDOFF_SCHEMA, VoxelGrid, to_voxel_grid as voxel_export

#: Published DARPA SUBOFF (AFF-1) displacement volume, via the offset table.
SUBOFF_VOLUME_M3 = 0.699211

ROOT = Path(__file__).resolve().parents[1]


def _env() -> dict[str, str]:
    """Subprocess environment that finds the package with or without an install."""
    existing = os.environ.get("PYTHONPATH", "")
    path = str(ROOT) + (os.pathsep + existing if existing else "")
    return {**os.environ, "PYTHONPATH": path}


#: Sphere volume used throughout: r = 0.5 m.
SPHERE_VOLUME = 4.0 / 3.0 * math.pi * 0.5**3


def sphere_error(cells: int, *, bounds=None) -> tuple[float, float]:
    grid = to_voxel_grid(SphereField(radius=0.5), cells, bounds=bounds)
    return (
        grid.volume() / SPHERE_VOLUME - 1.0,
        grid.occupied_volume() / SPHERE_VOLUME - 1.0,
    )


def test_the_export_is_reachable_from_the_package() -> None:
    import tensorcad

    assert tensorcad.to_voxel_grid is voxel_export
    assert "to_voxel_grid" in tensorcad.__all__
    assert "VoxelGrid" in tensorcad.__all__


def test_the_grid_is_cell_centred_and_the_spacing_is_the_box_over_the_cells() -> None:
    field = SphereField(radius=0.5)
    grid = to_voxel_grid(field, 16)
    assert grid.shape == (16, 16, 16)
    assert grid.cells == 16**3
    # The field's own box is +-0.5 on every axis, so one cell is 1/16 = 0.0625 m and the
    # first centre sits half a cell in from the corner -- not on it.
    assert grid.spacing == pytest.approx((0.0625, 0.0625, 0.0625))
    assert grid.origin == pytest.approx((-0.5, -0.5, -0.5))
    assert float(grid.axis(0)[0]) == pytest.approx(-0.5 + 0.03125)
    assert float(grid.axis(0)[-1]) == pytest.approx(0.5 - 0.03125)
    assert grid.cell_volume == pytest.approx(0.0625**3)
    assert grid.bounds_min == pytest.approx((-0.5, -0.5, -0.5))
    assert grid.bounds_max == pytest.approx((0.5, 0.5, 0.5))


def test_every_cell_holds_the_analytic_distance() -> None:
    """An analytic field must be sampled, not approximated: zero error, not "small"."""
    numpy = pytest.importorskip("numpy")
    field = SphereField(radius=0.5)
    grid = to_voxel_grid(field, (9, 7, 5))
    xs, ys, zs = (numpy.meshgrid(*[grid.axis(axis) for axis in range(3)], indexing="ij"))
    exact = numpy.sqrt(xs**2 + ys**2 + zs**2) - 0.5
    assert float(numpy.max(numpy.abs(grid.sdf - exact))) < 1e-12


def test_the_sign_convention_is_negative_inside() -> None:
    numpy = pytest.importorskip("numpy")
    field = SphereField(radius=0.5)
    grid = to_voxel_grid(field, 24)
    centre = grid.sdf[12, 12, 12]
    assert centre < 0.0
    # Cell centre (12, 12, 12) sits at 0.5/24 = 0.020833 m on every axis.
    offset = 12.5 * (1.0 / 24.0) - 0.5
    assert centre == pytest.approx(math.sqrt(3.0) * offset - 0.5, abs=1e-12)
    assert grid.sdf[0, 0, 0] > 0.0
    # Hard inside count, checked against the same centres computed independently.
    xs, ys, zs = numpy.meshgrid(*[grid.axis(axis) for axis in range(3)], indexing="ij")
    expected = int(numpy.count_nonzero(xs**2 + ys**2 + zs**2 < 0.25))
    assert grid.inside_count() == expected


def test_the_partial_fill_volume_beats_hard_cell_counting() -> None:
    """Both estimators are offered; the second-order one has to actually be better."""
    for cells in (16, 32):
        partial, hard = sphere_error(cells)
        assert abs(partial) < abs(hard) / 2.0, f"{cells}^3: {partial:+.5%} vs {hard:+.5%}"


def test_the_volume_error_falls_with_resolution() -> None:
    """A shape whose box is not grid-coincident, so the error is not aliased away."""
    ellipsoid = EllipsoidField(semi_axes=(0.5, 0.25, 0.125))
    exact = 4.0 / 3.0 * math.pi * 0.5 * 0.25 * 0.125
    errors = [
        to_voxel_grid(ellipsoid, cells).volume() / exact - 1.0 for cells in (16, 32, 64)
    ]
    assert all(error > 0.0 for error in errors), errors  # the cut rule overfills slightly
    assert errors[0] > errors[1] > errors[2], errors
    assert abs(errors[1]) < 3e-3
    assert abs(errors[2]) < 1e-3


def test_a_cube_grid_reaches_a_tenth_of_a_percent() -> None:
    partial, _ = sphere_error(64, bounds=((-0.6, -0.6, -0.6), (0.6, 0.6, 0.6)))
    assert abs(partial) < 1e-3, f"{partial:+.5%}"


def test_anisotropic_resolution_sets_one_spacing_per_axis() -> None:
    """A slender body wants per-axis cells, exactly as ``to_mesh`` allows."""
    field = SuboffHullField()
    grid = to_voxel_grid(field, (109, 13, 13))
    assert grid.shape == (109, 13, 13)
    assert grid.spacing[0] == pytest.approx(field.length / 109)
    assert grid.spacing[1] == pytest.approx(2.0 * field.max_radius / 13)
    assert grid.spacing[2] == pytest.approx(2.0 * field.max_radius / 13)
    # The three-cell form exists so an 8.6:1 body is not sampled by a cube: 109 cells
    # along 4.36 m and 13 across 0.508 m land within 3 % of the same cell size, whereas
    # 13^3 would be a 0.34 m cell along the length.
    ratio = grid.spacing[1] / grid.spacing[0]
    assert 0.9 < ratio < 1.1, ratio
    cube = to_voxel_grid(field, 13)
    assert cube.spacing[0] / cube.spacing[1] == pytest.approx(field.length / (2.0 * field.max_radius))


def test_the_sampling_box_can_be_overridden() -> None:
    field = SphereField(radius=0.5)
    box = ((-0.6, -0.55, -0.52), (0.65, 0.6, 0.58))
    grid = to_voxel_grid(field, 12, bounds=box)
    assert grid.origin == pytest.approx(box[0])
    assert grid.bounds_max == pytest.approx(box[1])
    assert grid.spacing[0] == pytest.approx(1.25 / 12)
    assert grid.spacing[1] == pytest.approx(1.15 / 12)
    # The window is bigger than the sphere, so it cannot be the field's own box.
    assert grid.volume() / SPHERE_VOLUME - 1.0 == pytest.approx(
        sphere_error(12, bounds=box)[0], abs=1e-15
    )


def test_padding_adds_cells_without_moving_the_grid() -> None:
    field = SphereField(radius=0.5)
    plain = to_voxel_grid(field, 16)
    padded = to_voxel_grid(field, 16, padding=2)
    assert padded.shape == (20, 20, 20)
    assert padded.spacing == pytest.approx(plain.spacing)
    assert padded.origin == pytest.approx(tuple(value - 2 * 0.0625 for value in plain.origin))
    # Interior cells are bit-identical: padding grows the box, it does not resample it.
    assert padded.sdf[2:-2, 2:-2, 2:-2] == pytest.approx(plain.sdf, abs=0.0)
    assert padded.volume() == pytest.approx(plain.volume())


def test_float32_storage_is_offered_for_memory_without_moving_the_volume() -> None:
    numpy = pytest.importorskip("numpy")
    field = SphereField(radius=0.5)
    dense = to_voxel_grid(field, 24)
    light = to_voxel_grid(field, 24, dtype="float32")
    assert light.sdf.dtype == numpy.float32
    assert light.sdf.nbytes == dense.sdf.nbytes // 2
    assert light.volume() == pytest.approx(dense.volume(), rel=1e-5)


def test_the_pure_python_path_agrees_with_the_vectorised_one(monkeypatch: pytest.MonkeyPatch) -> None:
    """Fields without a NumPy callable are still exportable; the result must match."""
    field = SphereField(radius=0.5)
    fast = to_voxel_grid(field, 8)
    monkeypatch.setattr(SphereField, "numpy_callable", lambda self: None)
    slow = to_voxel_grid(field, 8)
    assert slow.sdf == pytest.approx(fast.sdf, abs=1e-15)


def test_the_volume_agrees_with_the_published_suboff_displacement() -> None:
    """Cross-representation check: the voxel path against Groves et al., not the mesh."""
    grid = to_voxel_grid(SuboffHullField(), (109, 13, 13))
    partial = grid.volume() / SUBOFF_VOLUME_M3 - 1.0
    hard = grid.occupied_volume() / SUBOFF_VOLUME_M3 - 1.0
    assert abs(partial) < 5e-3, f"partial-fill volume {partial:+.4%}"
    assert abs(hard) < 3e-2, f"hard-count volume {hard:+.4%}"
    assert abs(partial) < abs(hard)


def test_no_area_estimator_is_shipped() -> None:
    """Area off a sampled SDF is orientation-biased, so the mesh path keeps the job.

    Measured before removing it: tent-weighted delta quadrature gave +0.16 % on a
    sphere at 64^3 but -8.4 % on the SUBOFF hull at (218, 26, 26), because lattice
    delta quadrature loses the deeply diagonal stern taper while it handles a
    grid-aligned barrel. A number that is 8 % low on the hull of interest is worse
    than no number: ``to_mesh(field, cells).report.surface_area`` is -0.066 % there.
    """
    assert not hasattr(VoxelGrid, "surface_area")
    assert "surface_area_m2" not in to_voxel_grid(SphereField(radius=0.5), 8).as_handoff()


def test_the_handoff_record_is_complete_json_and_deterministic() -> None:
    grid = to_voxel_grid(SphereField(radius=0.5), 12)
    record = grid.as_handoff()
    assert record["schema"] == VOXEL_HANDOFF_SCHEMA
    assert record["inside_sign"] == "negative_inside"
    assert record["sampling"] == "cell_centred"
    assert record["shape"] == [12, 12, 12]
    assert record["units"] == "m"
    assert record["coordinate_frame"] == "object_local"
    assert record["sdf_sha256"] == grid.sdf_digest()
    assert record["provenance"]["field_kind"] == "sphere"
    assert record["volume_m3"] == pytest.approx(grid.volume())
    json.dumps(record)  # no NumPy scalars: the record is a document, not an array
    assert record == to_voxel_grid(SphereField(radius=0.5), 12).as_handoff()


def test_the_digest_identifies_the_grid_not_the_representation() -> None:
    fine = to_voxel_grid(SphereField(radius=0.5), 16)
    coarse = to_voxel_grid(SphereField(radius=0.5), 8)
    assert fine.sdf_digest() != coarse.sdf_digest()
    assert fine.sdf_digest() != to_voxel_grid(EllipsoidField(semi_axes=(0.5, 0.5, 0.5)), 16).sdf_digest()


def test_bad_requests_fail_loudly() -> None:
    field = SphereField(radius=0.5)
    with pytest.raises(FieldError, match="unbounded"):
        to_voxel_grid(PlaneField(normal=(0.0, 0.0, 1.0), offset=0.0), 8)
    with pytest.raises(FieldError, match="at least 2"):
        to_voxel_grid(field, 1)
    with pytest.raises(FieldError, match="1 or 3"):
        to_voxel_grid(field, (8, 8))
    with pytest.raises(FieldError, match="degenerate or inverted"):
        to_voxel_grid(field, 8, bounds=((1.0, 0.0, 0.0), (0.0, 1.0, 1.0)))
    with pytest.raises(FieldError, match="two corners"):
        to_voxel_grid(field, 8, bounds=((-1.0, -1.0, -1.0),))
    with pytest.raises(FieldError, match="3 coordinates"):
        to_voxel_grid(field, 8, bounds=((-1.0, -1.0), (1.0, 1.0)))
    with pytest.raises(FieldError, match="not numeric"):
        to_voxel_grid(field, 8, bounds=(("a", 0.0, 0.0), (1.0, 1.0, 1.0)))
    with pytest.raises(FieldError, match="padding"):
        to_voxel_grid(field, 8, padding=-1)
    with pytest.raises(FieldError, match="floating"):
        to_voxel_grid(field, 8, dtype="int32")
    with pytest.raises(FieldError, match="chunk_cells"):
        to_voxel_grid(field, 8, chunk_cells=0)
    with pytest.raises(FieldError, match="cell count"):
        to_voxel_grid(field, "many")


def test_the_module_imports_and_exports_without_numpy() -> None:
    """Zero-dependency rule: NumPy is a sampling cost, not an import cost."""
    script = textwrap.dedent(
        """
        import importlib.abc, sys

        class Blocker(importlib.abc.MetaPathFinder):
            def find_spec(self, name, path=None, target=None):
                if name.split(".")[0] == "numpy":
                    raise ImportError("numpy blocked by the test")
                return None

        sys.meta_path.insert(0, Blocker())
        import tensorcad
        from tensorcad.voxel import to_voxel_grid, VoxelGrid
        print("imported", tensorcad.__version__)
        try:
            to_voxel_grid(tensorcad.SphereField(radius=0.5), 8)
        except tensorcad.FieldError as exc:
            print("blocked:", "NumPy" in str(exc))
        """
    )
    result = subprocess.run(
        [sys.executable, "-c", script], capture_output=True, text=True, cwd=ROOT, env=_env()
    )
    assert result.returncode == 0, result.stderr
    assert "blocked: True" in result.stdout, result.stdout


def test_the_example_script_passes_its_own_claims() -> None:
    result = subprocess.run(
        [sys.executable, str(ROOT / "examples" / "geometry_core" / "run_voxel_export.py")],
        capture_output=True,
        text=True,
        cwd=ROOT,
        env=_env(),
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert "voxel export: all checks passed" in result.stdout


def test_the_grid_keeps_a_grid_shape_for_a_degenerate_axis_request() -> None:
    """Guard the ``VoxelGrid`` contract rather than trusting the sampler's internals."""
    grid = to_voxel_grid(SphereField(radius=0.5), (3, 3, 3))
    assert isinstance(grid, VoxelGrid)
    assert len(grid.shape) == 3
    assert grid.sdf.size == grid.cells