"""Tests for the parametric hull model: topology, fidelity, validation, API."""
from __future__ import annotations

import dataclasses
import itertools
import json
import math

import pytest

from tensorcad import HullParameters, build_mesh, export_hull, generate_hull, validate_hull
from tensorcad import mesh as meshlib
from tensorcad.hull import BREADTH_EXPONENT, MAX_LONGITUDINAL_STATIONS, TAPER_INTEGRAL


def expected_vertex_count(nx: int, nz: int) -> int:
    """Interior stations carry a full section ring; the two ends collapse to a line."""
    return (nx - 2) * (2 * nz - 1) + 2 * nz


def expected_face_count(nx: int, nz: int) -> int:
    """Lofted shell plus the waterline cap, minus the two collapsed cap quads."""
    return 4 * (nx - 1) * (nz - 1) + 2 * nx - 4


def analytic_volume(length: float, beam: float, draft: float) -> float:
    return beam * length * draft * TAPER_INTEGRAL * (2.0 / 3.0)


@pytest.mark.parametrize("nx,nz", [(3, 3), (5, 3), (11, 5), (41, 13), (12, 7)])
def test_topology_counts_are_exact(nx: int, nz: int) -> None:
    model = generate_hull(HullParameters(nx=nx, nz=nz))
    assert len(model.vertices) == expected_vertex_count(nx, nz)
    assert len(model.faces) == expected_face_count(nx, nz)
    assert len(model.offsets) == nx


@pytest.mark.parametrize("nx,nz", [(3, 3), (5, 4), (11, 5), (41, 13), (60, 9)])
def test_generated_surface_is_a_closed_oriented_solid(nx: int, nz: int) -> None:
    """The regression this package exists to fix: a mesh a solver can accept."""
    model = generate_hull(HullParameters(nx=nx, nz=nz))
    assert meshlib.is_watertight(model.mesh), "surface must be closed"
    assert not meshlib.boundary_edges(model.mesh)
    assert meshlib.is_edge_manifold(model.mesh)
    assert not meshlib.zero_area_faces(model.mesh)
    assert meshlib.total_signed_volume(model.mesh) > 0, "faces must point outward"
    edges = len(meshlib.edge_face_counts(model.mesh))
    assert len(model.vertices) - edges + len(model.faces) == 2, "closed genus-0 shell"


def test_default_hull_is_deterministic() -> None:
    first = generate_hull()
    second = generate_hull()
    assert first == second
    assert validate_hull(first) == []


def test_waterline_lies_exactly_on_z_zero() -> None:
    """Non-cancelling floats here would silently break coplanarity checks."""
    model = generate_hull(HullParameters(draft=0.8, nz=13))
    assert model.mesh.bounding_box().maximum[2] == 0.0
    assert model.mesh.bounding_box().minimum[2] == pytest.approx(-0.8, abs=0.0)
    keel = [point for point in model.vertices if point[2] == -0.8]
    assert keel and all(point[1] == 0.0 for point in keel), "the keel is one shared line"


def test_dimensions_change_geometry() -> None:
    baseline = generate_hull(HullParameters(length=10, beam=2, draft=0.8, nx=11, nz=5))
    wider = generate_hull(HullParameters(length=10, beam=3, draft=0.8, nx=11, nz=5))
    deeper = generate_hull(HullParameters(length=10, beam=2, draft=1.2, nx=11, nz=5))
    assert max(y for _, y, _ in wider.vertices) == pytest.approx(1.5)
    assert max(y for _, y, _ in wider.vertices) > max(y for _, y, _ in baseline.vertices)
    assert deeper.mesh.bounding_box().minimum[2] == pytest.approx(-1.2)
    assert wider.hydrostatics.displacement_volume > baseline.hydrostatics.displacement_volume
    assert deeper.hydrostatics.displacement_volume > baseline.hydrostatics.displacement_volume


def test_offsets_describe_the_waterline_and_close_at_both_ends() -> None:
    model = generate_hull(HullParameters(beam=2.4, nx=41, nz=13))
    assert model.offsets[0][1] == 0.0
    assert model.offsets[-1] == (10.0, 0.0)
    assert max(width for _, width in model.offsets) == pytest.approx(2.4)
    # For this form the widest point of every section is on the waterline.
    for x, width in model.offsets:
        row = max(abs(y) for px, y, _ in model.vertices if math.isclose(px, x, abs_tol=1e-9))
        assert width == pytest.approx(2.0 * row)


def test_even_station_count_reports_the_beam_shortfall_instead_of_hiding_it() -> None:
    odd = generate_hull(HullParameters(nx=41, nz=13))
    even = generate_hull(HullParameters(nx=40, nz=13))
    assert odd.validation.warnings == []
    assert any("not_achievable" in warning for warning in even.validation.warnings)
    realised = max(width for _, width in even.offsets)
    assert realised < 2.0
    assert even.validation.passed, "a sampling deficit is a warning, not a broken model"


def test_hydrostatics_match_the_analytic_surface_within_grid_error() -> None:
    for nx, nz, tolerance in [(11, 5, 0.06), (41, 13, 0.015), (161, 41, 0.002)]:
        model = generate_hull(HullParameters(nx=nx, nz=nz))
        exact = analytic_volume(10.0, 2.0, 0.8)
        assert model.hydrostatics.displacement_volume == pytest.approx(exact, rel=tolerance)


def test_hydrostatics_are_consistent_with_each_other() -> None:
    model = generate_hull(HullParameters(length=12, beam=2.4, draft=0.9, nx=41, nz=13))
    hydro = model.hydrostatics
    box = model.mesh.bounding_box()
    assert hydro.block_coefficient == pytest.approx(
        hydro.displacement_volume / (box.sizes[0] * box.sizes[1] * box.sizes[2])
    )
    assert hydro.displacement_mass == pytest.approx(hydro.displacement_volume * 1025.0)
    # A fore-aft symmetric form must have its buoyancy and flotation centres aligned.
    assert hydro.longitudinal_centre_buoyancy == pytest.approx(hydro.longitudinal_centre_flotation, rel=1e-9)
    assert hydro.longitudinal_centre_buoyancy == pytest.approx(6.0, rel=1e-9)
    assert hydro.vertical_centre_buoyancy < 0.0
    assert hydro.waterplane_area > 0.0 and hydro.wetted_surface_area > hydro.waterplane_area


def test_resolution_refinement_reduces_the_error() -> None:
    errors = []
    for nx, nz in [(11, 5), (21, 9), (41, 13), (81, 21)]:
        model = generate_hull(HullParameters(nx=nx, nz=nz))
        errors.append(abs(model.hydrostatics.displacement_volume / analytic_volume(10, 2, 0.8) - 1.0))
    assert all(a > b for a, b in itertools.pairwise(errors))


def test_discretisation_estimate_reports_measured_change() -> None:
    from tensorcad import estimate_discretisation_error

    report = estimate_discretisation_error(HullParameters(nx=21, nz=9))
    assert report["status"] == "measured"
    assert report["refined_grid"] == {"nx": 41, "nz": 17}
    for value in report["relative_change"].values():
        assert 0.0 < value < 0.1
    skipped = estimate_discretisation_error(HullParameters(nx=1500, nz=270))
    assert skipped["status"] == "skipped"


@pytest.mark.parametrize(
    "kwargs",
    [
        {"length": 0},
        {"beam": -1},
        {"draft": float("nan")},
        {"draft": float("inf")},
        {"nx": 2},
        {"nz": 1},
        {"nx": 11.0},
        {"nz": True},
        {"nx": 10**8},
        {"length": 10**9},
    ],
)
def test_invalid_parameters_are_rejected(kwargs: dict) -> None:
    values: dict = {"length": 10.0, "beam": 2.0, "draft": 0.8, "nx": 11, "nz": 5}
    values.update(kwargs)
    with pytest.raises(ValueError):
        generate_hull(HullParameters(**values))  # type: ignore[arg-type]


def test_all_parameter_problems_are_reported_at_once() -> None:
    errors, _ = HullParameters(length=-1, beam=0, draft=float("nan"), nx=2, nz=2).problems()
    assert len(errors) == 5, "one message per bad value, not first-failure-only"


def test_implausible_dimension_ratios_warn_without_failing() -> None:
    model = generate_hull(HullParameters(length=1.0, beam=50.0, draft=0.8, nx=7, nz=5))
    assert model.validation.passed
    assert any("length_to_beam" in warning for warning in model.validation.warnings)
    errors, warnings = HullParameters(length=1.0, beam=50.0).problems()
    assert not errors and warnings


def test_station_grid_limits_are_enforced_before_allocating() -> None:
    assert MAX_LONGITUDINAL_STATIONS == 4001
    with pytest.raises(ValueError, match="nx must be <="):
        generate_hull(HullParameters(nx=MAX_LONGITUDINAL_STATIONS + 1))


def test_grid_symmetry_flag_matches_the_station_parity() -> None:
    assert HullParameters(nx=41).grid_is_symmetric
    assert not HullParameters(nx=40).grid_is_symmetric


def test_generation_refuses_to_hand_back_an_unusable_mesh(monkeypatch) -> None:
    """A broken surface must raise, never be exported with a clean bill of health."""
    from tensorcad import hull as hull_module
    from tensorcad.geometry import Mesh

    original = hull_module.build_mesh

    def broken(parameters: HullParameters) -> tuple[Mesh, tuple[tuple[float, float], ...]]:
        good, offsets = original(parameters)
        # Punch a hole in an otherwise valid mesh: closed -> open.
        return Mesh(vertices=good.vertices, faces=good.faces[:-2]), offsets

    monkeypatch.setattr(hull_module, "build_mesh", broken)
    with pytest.raises(ValueError, match="generated hull failed validation"):
        hull_module.generate_hull(HullParameters(nx=11, nz=5))


def test_mesh_is_built_from_a_single_source_of_truth() -> None:
    """build_mesh is what generate_hull uses, so both must agree exactly."""
    parameters = HullParameters(length=8.0, beam=1.6, draft=0.7, nx=17, nz=7)
    mesh, offsets = build_mesh(parameters)
    model = generate_hull(parameters)
    assert model.mesh.vertices == mesh.vertices
    assert model.mesh.faces == mesh.faces
    assert model.offsets == offsets


def test_metadata_is_complete_and_json_serialisable(tmp_path) -> None:
    model = generate_hull(HullParameters(nx=11, nz=5))
    metadata = model.metadata(with_discretisation=True)
    text = json.dumps(metadata, sort_keys=True, allow_nan=False)
    loaded = json.loads(text)
    assert loaded["schema_version"] == "2.0"
    assert loaded["units"] == "m"
    assert "coordinate_frame" in loaded
    assert loaded["parameters"] == {
        "length": 10.0, "beam": 2.0, "draft": 0.8, "nx": 11, "nz": 5, "water_density": 1025.0,
    }
    assert loaded["validation"]["passed"] is True
    assert loaded["mesh"]["vertex_count"] == len(model.vertices)
    assert loaded["hydrostatics"]["displacement_volume_m3"] > 0
    assert loaded["discretisation"]["status"] == "measured"
    assert "not implemented" in loaded["engineering_status"]
    assert loaded["generator"]["longitudinal_taper"] == f"sin(pi*x/length)^{BREADTH_EXPONENT}"


def test_export_writes_reproducible_round_trippable_artifacts(tmp_path) -> None:
    from tensorcad import read_metadata, read_obj

    model = generate_hull(HullParameters(nx=9, nz=5))
    artifacts = export_hull(model, tmp_path, convergence=True)
    assert set(artifacts) == {"mesh", "offsets", "metadata"}
    assert all((tmp_path / name).is_file() for name in ("hull.obj", "offsets.csv", "metadata.json"))

    written = read_metadata(artifacts["metadata"])
    assert written["validation"]["errors"] == []
    assert written["units"] == "m"
    assert written["mesh"]["vertex_count"] == 7 * 9 + 10

    # Exact round-trip: the OBJ must reload as the same model, not an approximations of it.
    reloaded = read_obj(artifacts["mesh"])
    assert reloaded.vertices == model.vertices
    assert reloaded.faces == model.faces
    assert meshlib.is_watertight(reloaded)

    rows = (tmp_path / "offsets.csv").read_text().splitlines()
    assert rows[0] == "x_m,overall_breadth_m"
    assert len(rows) == 1 + 9
    assert next(float(row.split(",")[1]) for row in rows[1:]) == 0.0


def test_parameters_are_immutable_and_positional_compatible() -> None:
    legacy = HullParameters(10.0, 2.0, 0.8, 41, 13)  # the 0.1 positional form still works
    assert legacy.water_density == 1025.0
    assert legacy == HullParameters()
    with pytest.raises(dataclasses.FrozenInstanceError):  # inputs cannot change under a model
        legacy.length = 1.0  # type: ignore[misc]
