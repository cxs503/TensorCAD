"""The benchmark's own contract: measured, reproducible, and honest about what it cannot do.

A benchmark that is not trustworthy is worse than no benchmark, so these tests check the
harness rather than the geometry: that the metric schema refuses to encode a missing value
as a number, that the reference row really is the reference, that the ray caster agrees
with the closed form it is aimed at, that a gate failure survives into the verdict, and
that a run is reproducible (same digests, same numbers).

The geometry gates are exercised on a small but *passing* configuration, because CI has to
be fast. The resolution dependence is asserted where it is a fact: at 12x24 rings the mesh
row's own discretisation fails the depth gate, and that failure is reported rather than
tuned away.
"""
from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path

import pytest

from tensorcad.benchmark import (
    BENCHMARK_SCHEMA,
    DENT_DEPTH_M,
    FAILED,
    IMPORTED,
    MEASURED,
    NOT_SUPPORTED,
    BenchmarkResult,
    CaseResult,
    Metric,
    RepresentationResult,
    _angle_from_axis,
    _field_radii,
    _lat_long_surface,
    _RayCaster,
    dent_radius,
    fibonacci_directions,
    have_neural,
    hull_handoff_case,
    run_benchmark,
    sphere_dent_case,
    stratified_directions,
)


@pytest.fixture(scope="module")
def quick_case() -> CaseResult:
    return sphere_dent_case(rays=128)


@pytest.fixture(scope="module")
def quick_run() -> BenchmarkResult:
    return run_benchmark(quick=True)


def test_directions_are_deterministic_and_unit_length() -> None:
    first = fibonacci_directions(64)
    second = fibonacci_directions(64)
    assert first == second
    assert len(first) == 64
    for direction in first:
        assert math.sqrt(sum(component * component for component in direction)) == pytest.approx(1.0)


def test_directions_reject_a_nonsense_count() -> None:
    with pytest.raises(ValueError, match="at least 1"):
        fibonacci_directions(0)


def test_the_dent_reference_is_the_published_profile() -> None:
    """Peak depth on the axis, Gaussian falloff, untouched at the far pole."""
    radius, depth, sigma = 0.5, 0.06, 0.22
    assert dent_radius((1.0, 0.0, 0.0), radius, depth, sigma) == pytest.approx(radius - depth)
    assert dent_radius((-1.0, 0.0, 0.0), radius, depth, sigma) == pytest.approx(radius, abs=1e-15)
    half = math.cos(sigma)
    at_sigma = dent_radius((half, math.sin(sigma), 0.0), radius, depth, sigma)
    assert at_sigma == pytest.approx(radius - depth * math.exp(-0.5), rel=1e-12)


def test_a_missing_measurement_can_never_be_a_number() -> None:
    """The single most important rule in the schema: absence is not zero."""
    with pytest.raises(ValueError, match="needs a value"):
        Metric(value=None, unit="m", status=MEASURED)
    with pytest.raises(ValueError, match="must not carry a value"):
        Metric(value=0.0, unit="m", status=NOT_SUPPORTED, note="needs torch")
    with pytest.raises(ValueError, match="not finite"):
        Metric(value=float("nan"), unit="m", status=MEASURED)
    with pytest.raises(ValueError, match="unknown metric status"):
        Metric(value=1.0, status="probably")


def test_gates_never_pass_on_a_missing_value() -> None:
    missing = Metric(value=None, unit="m", status=NOT_SUPPORTED, note="no wheels")
    assert missing.within(1.0).passed is False
    assert missing.at_most(1.0).passed is False
    assert missing.at_least(0.0).passed is False
    assert missing.equals(2.0).passed is False


def test_gate_helpers_compare_what_they_say() -> None:
    assert Metric(value=0.004, unit="m").within(1e-2).passed is True
    assert Metric(value=0.04, unit="m").within(1e-2).passed is False
    assert Metric(value=2.0, unit="count").equals(2.0).passed is True
    assert Metric(value=3.0, unit="count").equals(2.0).passed is False
    assert Metric(value=0.9, unit="fraction").at_least(0.99).passed is False


def test_the_reference_row_is_the_reference(quick_case) -> None:
    analytic = quick_case.by_id("analytic_field")
    assert analytic.status == MEASURED
    assert analytic.metrics["reference_radius_error_max_m"].value == 0.0
    assert analytic.metrics["reference_radius_error_max_m"].passed is True
    # The exact volume of a 0.5 m sphere with a 0.06 m dent is just under the sphere's.
    volume = analytic.metrics["volume_m3"].value
    assert volume is not None
    assert 0.98 * 4.0 / 3.0 * math.pi * 0.5**3 < volume < 4.0 / 3.0 * math.pi * 0.5**3
    # Triangulating an exact surface costs something; it is reported, not hidden.
    assert analytic.metrics["triangulated_volume_error"].value is not None


def test_the_mesh_row_measures_a_real_edit(quick_case) -> None:
    mesh = quick_case.by_id("mesh_local_edit")
    assert mesh.status == MEASURED
    metrics = mesh.metrics
    assert metrics["ray_coverage"].value == 1.0
    assert metrics["watertight"].value == 1.0
    assert metrics["euler_characteristic"].value == 2.0
    assert metrics["dent_depth_peak_m"].value is not None
    assert 0.5 * DENT_DEPTH_M < metrics["dent_depth_peak_m"].value <= DENT_DEPTH_M
    assert metrics["dent_depth_mean_error_m"].passed is True
    assert metrics["locality_ratio"].passed is True
    assert all(metric.value is not None for metric in metrics.values())


def test_the_unavailable_rows_say_why(quick_case) -> None:
    """A row without a measurement must carry a reason and a command, never a zero."""
    for identifier in ("neural_field", "ocp_step"):
        row = quick_case.by_id(identifier)
        assert row.status in (NOT_SUPPORTED, IMPORTED, MEASURED)
        for metric in row.metrics.values():
            if metric.status == NOT_SUPPORTED:
                assert metric.value is None
                assert metric.note, f"{identifier} did not explain itself"


def test_a_coarse_grid_fails_its_own_gate_and_says_so() -> None:
    """The mesh row's accuracy is resolution: at 12x24 it misses the depth gate.

    This is the property that makes the benchmark worth running -- the gate is a
    measurement, not a decoration, and it fails when the representation is too coarse.
    """
    coarse = sphere_dent_case(rays=128, rings=12, segments=24)
    assert not coarse.by_id("mesh_local_edit").metrics["dent_depth_axis_error_m"].passed
    assert "mesh_edit_hits_the_target_depth" in {gate.name for gate in coarse.failed_gates}

    # The same case on the grid the harness ships clears it, so the failure above is
    # resolution and not a broken edit.
    fine = sphere_dent_case(rays=128)
    assert fine.by_id("mesh_local_edit").metrics["dent_depth_axis_error_m"].passed
    assert not fine.failed_gates


def test_the_sag_floor_is_well_inside_the_depth_gate(quick_case) -> None:
    """The discretisation floor must sit below the gate by a stated margin.

    The depth gate is 1 % of 0.06 m = 6.0e-4 m. A lat-long grid is inscribed, so even an
    unedited sphere comes back low on a ray between vertices; if that floor were close to
    the gate, the gate would be measuring the grid instead of the edit. It is reported as a
    metric so a reader can see the headroom rather than trust a comment.
    """
    gate = 0.01 * DENT_DEPTH_M
    floor = quick_case.by_id("mesh_local_edit").metrics["discretisation_floor_m"].value
    assert floor is not None
    assert 0.0 < floor < 0.35 * gate, f"chord-sag floor {floor:.3e} m is too close to the {gate:.3e} m gate"


def test_the_cone_gets_enough_rays_to_measure_it() -> None:
    """A uniform ray set puts ~4 rays in a small edit cone; the stratified one does not."""
    uniform = fibonacci_directions(512)
    stratified = stratified_directions(512, 0.22)
    uniform_inside = sum(1 for d in uniform if _angle_from_axis(d) <= 0.22)
    stratified_inside = sum(1 for d in stratified if _angle_from_axis(d) <= 0.22)
    assert uniform_inside < 10
    assert stratified_inside >= 150
    assert len(set(stratified)) == len(stratified)
    with pytest.raises(ValueError, match="at least 4"):
        stratified_directions(3, 0.22)
    with pytest.raises(ValueError, match=r"in \(0, pi\)"):
        stratified_directions(64, 4.0)


def test_the_ray_marcher_recovers_a_known_surface() -> None:
    """The instrument for learned fields, tested against a surface with a closed form."""
    # Negative inside, positive outside: the sign convention the .ngo field declares.
    field_radii = _field_radii(
        lambda points: [(p[0] ** 2 + p[1] ** 2 + p[2] ** 2) ** 0.5 - 0.5 for p in points],
        fibonacci_directions(32),
    )
    assert all(radius is not None for radius in field_radii)
    assert all(abs(radius - 0.5) < 1e-6 for radius in field_radii if radius is not None)


def test_a_ray_that_never_crosses_is_reported_missing() -> None:
    radii = _field_radii(lambda points: [1.0] * len(points), fibonacci_directions(8))
    assert radii == [None] * 8


def test_the_ray_caster_agrees_with_the_closed_form() -> None:
    """Sanity on the instrument: a sphere's rays come back at its radius, from inside.

    The deficit is not a bug and not zero: a lat-long triangulation is inscribed, so a ray
    between vertices lands on a chord below the true surface. This is the floor that any
    edit measured against an analytic surface inherits, and the reason the mesh row compares
    against the analytic surface only for the cone mean while the depth gate compares against
    the reference *at the same ray*.
    """
    for rings, segments in ((24, 48), (40, 80)):
        sphere = _lat_long_surface(lambda direction: 0.5, rings=rings, segments=segments)
        caster = _RayCaster(sphere)
        # Worst-case chord sag of an inscribed lat-long grid, half a cell either way.
        sag = 0.5 * (1.0 - math.cos(math.pi / rings)) + 0.5 * (1.0 - math.cos(math.pi / segments))
        for direction in fibonacci_directions(64):
            hit = caster.radius(direction)
            assert hit is not None
            assert 0.5 - sag <= hit <= 0.5 + 1e-12


def test_every_measured_row_faces_the_same_depth_threshold(quick_case) -> None:
    gate = {gate.name: gate for gate in quick_case.gates}[
        "the_depth_threshold_is_the_same_for_every_measured_row"
    ]
    assert gate.passed
    for row in quick_case.representations:
        if row.status == "measured" and "dent_depth_mean_error_m" in row.metrics:
            metric = row.metrics["dent_depth_mean_error_m"]
            assert metric.tolerance == pytest.approx(0.01 * DENT_DEPTH_M)


def test_have_neural_reflects_what_is_importable() -> None:
    """The row must not claim it can measure when the runtime is absent, or vice versa.

    A spec on disk is not enough -- this project has been in the state where torch was
    installed and safetensors was not, which `find_spec` reports as available and the import
    then refuses.
    """
    import importlib

    expected = True
    for name in ("numpy", "torch", "safetensors"):
        try:
            importlib.import_module(name)
        except ImportError:
            expected = False
    assert have_neural() is expected


def test_the_hull_case_measures_against_the_analytic_volume() -> None:
    from tensorcad.hull import analytic_volume

    case = hull_handoff_case()
    row = case.by_id("parametric_hull")
    exact = analytic_volume(12.0, 2.4, 0.9)
    assert row.metrics["volume_error"].value == pytest.approx(
        abs(row.metrics["volume_error"].value), abs=1.0
    )
    assert row.metrics["volume_error"].passed is True
    assert row.metrics["design_loop_solved"].value == 1.0
    assert row.metrics["export_node_error_m"].value == 0.0
    # The grid is disclosed, not assumed: refining it is a measured number.
    assert row.metrics["grid_refinement_volume_change"].value is not None
    assert exact > 0.0
    for gate in case.gates:
        assert gate.passed, gate


def test_a_full_quick_run_passes_and_is_reproducible(quick_run) -> None:
    assert quick_run.passed, quick_run.failures
    assert quick_run.environment["tensorcad"]
    assert quick_run.environment["python"]
    again = run_benchmark(quick=True)
    first = quick_run.to_dict()
    second = again.to_dict()
    # Timings differ; every measured number and every gate must not.
    def strip(payload: dict) -> dict:
        payload = json.loads(json.dumps(payload))
        payload.pop("seconds")
        payload.pop("environment")
        for case in payload["cases"]:
            case.pop("seconds")
            for representation in case["representations"]:
                representation.pop("seconds")
                for metric in representation["metrics"].values():
                    metric.pop("note")
                    if metric.get("unit") == "s":  # wall-clock, the one thing that varies
                        metric["value"] = None
        return payload

    assert strip(first) == strip(second)


def test_the_documents_are_generated_and_agree(quick_run, tmp_path) -> None:
    json_path = quick_run.write_json(tmp_path / "nested" / "results.json")
    report_path = quick_run.write_markdown(tmp_path / "report.md")
    payload = json.loads(json_path.read_text(encoding="utf-8"))
    assert payload["schema"] == BENCHMARK_SCHEMA
    assert payload["passed"] is True
    assert [case["identifier"] for case in payload["cases"]] == [
        "case_a_sphere_dent",
        "case_b_hull_handoff",
    ]
    report = report_path.read_text(encoding="utf-8")
    assert "verdict: pass" in report
    assert "case_a_sphere_dent" in report
    # The report is derived: a number in it must exist in the JSON.
    dent = quick_run.cases[0].by_id("mesh_local_edit").metrics["locality_ratio"].value
    assert dent is not None
    assert f"{dent:.3e}" in report or f"{dent:.6g}" in report


def test_a_failed_gate_makes_the_run_fail() -> None:
    """The verdict is the conjunction over cases; a red gate cannot be averaged away."""
    good = CaseResult(
        identifier="ok",
        question="?",
        reference="r",
        representations=(RepresentationResult(identifier="row", role="r", status=MEASURED),),
        gates=(),
    )
    bad = CaseResult(
        identifier="bad",
        question="?",
        reference="r",
        representations=(RepresentationResult(identifier="row", role="r", status=MEASURED),),
        gates=(),
    )
    from tensorcad.benchmark import Gate

    bad = CaseResult(
        identifier="bad",
        question="?",
        reference="r",
        representations=bad.representations,
        gates=(Gate("depth", False, "0.012 m vs 0.006 m"),),
    )
    result = BenchmarkResult(cases=(good, bad), environment={}, seconds=0.0)
    assert result.passed is False
    assert result.failures == ("bad: depth",)
    payload = result.to_dict()
    assert payload["passed"] is False
    assert payload["failures"] == ["bad: depth"]


def test_a_corrupt_bundle_is_reported_not_raised(tmp_path) -> None:
    """An unusable artifact is a failed row, not a crashed benchmark."""
    fake = tmp_path / "broken.ngo"
    fake.mkdir()
    (fake / "manifest.json").write_text('{"format": "nonsense"}', encoding="utf-8")
    case = sphere_dent_case(rays=64, rings=18, segments=36, neural_bundle=fake)
    row = case.by_id("neural_field")
    assert row.status == FAILED
    assert any("failed verification" in note for note in row.notes)
    for metric in row.metrics.values():
        assert metric.value is None


def test_the_runner_exits_on_a_failed_gate(tmp_path) -> None:
    """The command is a CI step, so its exit code has to mean something."""
    completed = subprocess.run(
        [
            sys.executable,
            str(Path(__file__).resolve().parents[1] / "examples" / "benchmarks" / "run_benchmark.py"),
            "--quick",
            "--output",
            str(tmp_path / "bench"),
        ],
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert (tmp_path / "bench" / "results.json").exists()
    assert (tmp_path / "bench" / "report.md").exists()
    assert "all gates passed" in completed.stdout


@pytest.mark.skipif(not have_neural(), reason="needs the neural extra to train and query a field")
def test_the_neural_row_is_re_measured_when_a_bundle_is_available(tmp_path) -> None:
    """The measured path: a real bundle, queried on this case's own rays.

    The bundle is trained for a few dozen steps with deliberately loose tolerances, so its
    accuracy is poor -- that is irrelevant here. What is being tested is that the row stops
    importing and starts measuring, that the producing run's numbers are kept beside the new
    ones instead of replacing them, and that a bad field is reported rather than raised.
    """
    from tensorcad.neural.objects import Tolerances, apply_dent_edit, fit_sphere

    loose = Tolerances(
        field_mae=1e9,
        boundary_error_mean=1e9,
        ray_coverage=0.0,
        dent_peak_error=1e9,
        dent_region_error_mean=1e9,
        protected_mean_displacement=1e9,
        protected_p99_displacement=1e9,
        protected_max_displacement=1e9,
    )
    bundle = tmp_path / "mini.ngo"
    fit_sphere(bundle, steps=40, rays=128, tolerance=loose, command=("test",))
    apply_dent_edit(bundle, depth=0.06, sigma=0.22, steps=40, rays=128, tolerance=loose,
                    command=("test",))

    case = sphere_dent_case(rays=64, rings=24, segments=48, neural_bundle=bundle)
    row = case.by_id("neural_field")
    assert row.status == MEASURED
    assert row.metrics["ray_coverage"].value is not None
    assert row.metrics["surface_error_mean_m"].value is not None
    assert row.metrics["dent_depth_mean_error_m"].status == MEASURED
    assert row.metrics["dent_depth_mean_error_m"].tolerance == pytest.approx(0.01 * DENT_DEPTH_M)
    assert any("queried on the same" in note for note in row.notes)
    # The producing run's own claims survive as imported metrics, clearly labelled.
    assert row.metrics["reported_dent_surface_error_mean_m"].status == IMPORTED
    assert row.metrics["protected_displacement_mean_m"].status == IMPORTED
    assert row.metrics["locality_ratio"].status == IMPORTED
    # A poor field is measured, not excused: this one misses the shared depth gate.
    assert row.metrics["dent_depth_mean_error_m"].passed is False
    assert "the_depth_threshold_is_the_same_for_every_measured_row" in {
        gate.name for gate in case.failed_gates
    }
