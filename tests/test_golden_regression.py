"""Tests for the golden geometry regression suite (``tensorcad.golden``).

Two layers are tested, on purpose:

* the **engine** — blessing, drift detection, hand-edited references, missing baselines,
  unknown cases, derived tolerances, fail-closed behaviour — using a synthetic case whose
  numbers are a formula, so the tests stay fast and deterministic;
* the **golden library** — the four committed cases (Wigley, Series 60, sphere/ellipsoid,
  DTMB 5415) checked against the committed baseline, which is the gate itself.

The synthetic case is deliberately *not* a hull: the engine must not depend on ship
semantics, and a test that reuses the real cases cannot tell an engine defect from a
change in the geometry pipeline.
"""

from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path

import pytest

from tensorcad import golden

ROOT = Path(__file__).resolve().parents[1]
BASELINE = ROOT / "docs" / "golden" / "tensorcad_geometry_golden.json"


class ProbeCase(golden.GoldenCase):
    """A case whose reference is a formula and whose measurement can be perturbed."""

    name = "probe"
    title = "synthetic probe (reference = 4/3·π·r³ with r = 1)"
    reference_source = "tests.test_golden_regression.ProbeCase.reference (closed form)"

    def __init__(self, scale: float = 1.0, invariant_residual: float = 0.0) -> None:
        self.scale = scale
        self.invariant_residual = invariant_residual

    def _truth(self) -> float:
        return 4.0 / 3.0 * math.pi * 1.0**3

    def measure(self, level: str) -> dict[str, float]:
        offset = {"coarse": 3e-3, "fine": 8e-4}[level]
        return {
            "value": self._truth() * (1.0 + offset) * self.scale,
            "coarse_only": float({"coarse": 1.0, "fine": 2.0}[level]),
        }

    def reference(self) -> dict[str, float]:
        return {"value": self._truth()}

    def invariants(self, values: dict[str, float]) -> dict[str, tuple[float, float, str]]:
        return {"identity_probe": (self.invariant_residual, 1e-9, "probe identity")}


def _bless(tmp_path: Path, case: golden.GoldenCase) -> Path:
    path = tmp_path / "baseline.json"
    golden.bless(path, cases=[case])
    return path


def _codes(result: golden.SuiteResult) -> set[str]:
    return {finding.code for finding in result.all_findings}


# --------------------------------------------------------------------------------------
# the committed golden library
# --------------------------------------------------------------------------------------
@pytest.fixture(scope="module")
def committed() -> golden.SuiteResult:
    return golden.run_suite()


def test_committed_baseline_passes(committed: golden.SuiteResult) -> None:
    """The gate the CI workflow and the git hook run: the tree must match its baseline."""
    assert committed.passed, [f"{f.case}/{f.code}: {f.detail}" for f in committed.errors]
    assert not committed.errors
    assert {case.name for case in committed.cases} == {
        "wigley",
        "series60",
        "sphere_ellipsoid",
        "dtmb5415",
    }


def test_golden_library_covers_the_required_forms(committed: golden.SuiteResult) -> None:
    """The four forms the anchor names: 5415_exact, Wigley, Series 60, sphere-ellipsoid."""
    cases = {case.name: case for case in committed.cases}
    assert "closed form" in cases["wigley"].reference_source
    assert "analytic_*" in cases["series60"].reference_source
    assert "π" in cases["sphere_ellipsoid"].reference_source
    assert "npz" in cases["dtmb5415"].reference_source
    for case in committed.cases:
        assert case.metrics, f"{case.name} asserts nothing"


def test_dtmb5415_publishes_topology_instead_of_a_fake_displacement(
    committed: golden.SuiteResult,
) -> None:
    """The 5415 shell is an open sheet: its tetrahedron sum must never be sold as a volume.

    Regression for the recorded "surface-integral volume disagrees with the cloud integral"
    limitation: the 2.6 relative discrepancy was an artifact of forcing ``closed=True`` on
    a mesh that is open at the deck edge, not an integration error.
    """
    case = {result.name: result for result in committed.cases}["dtmb5415"]
    asserted = {metric.name for metric in case.metrics}
    assert "surface_integral_volume_m3" not in asserted
    measured = case.measurements["fine"]
    assert {"mesh_is_watertight", "boundary_edge_count", "boundary_component_count",
            "signed_tetra_sum_open_shell_m3"} <= set(measured)
    assert measured["mesh_is_watertight"] == 0.0
    assert measured["boundary_edge_count"] > 0.0
    assert measured["boundary_component_count"] == 1.0
    assert "identity_open_shell_is_not_a_solid" in case.invariants
    limitations = {entry["id"]: entry for entry in case.known_limitations}
    assert "displacement-withheld-open-shell" in limitations
    withheld = limitations["displacement-withheld-open-shell"]
    assert withheld["asserted"] is False
    assert withheld["signed_tetra_sum_open_shell_m3"] < 0.0 < withheld["displacement_volume_m3"]


def test_forcing_closed_on_the_open_5415_shell_is_visible() -> None:
    """``closed=True`` on an open shell is a caller error; the report must expose it."""
    from tensorcad import mesh as meshlib
    from tensorcad.dtmb5415 import DTMB5415Space, WATERLINE_Z
    from tensorcad.hydrostatics import compute_hydrostatics

    space = DTMB5415Space()
    mesh = space.mesh()
    assert not meshlib.is_watertight(mesh)
    auto = compute_hydrostatics(mesh, waterline_z=WATERLINE_Z)
    assert "not trustworthy" in auto.basis
    forced = compute_hydrostatics(mesh, waterline_z=WATERLINE_Z, closed=True)
    assert forced.basis == "closed surface integrals"
    # The forced number is the signed tetrahedron sum of an open surface, not a volume:
    # it is negative while the displacement reported by the hull model is positive.
    assert forced.displacement_volume < 0.0 < space.hydrostatics()["displacement_volume_m3"]


def test_every_analytic_metric_carries_a_derived_tolerance(
    committed: golden.SuiteResult,
) -> None:
    analytic = [m for case in committed.cases for m in case.metrics if m.kind == "analytic"]
    assert analytic
    for metric in analytic:
        assert metric.tolerance is not None and metric.tolerance > 0.0
        assert "two-grid" in metric.basis or "floor" in metric.basis
        assert metric.reference is not None


def test_measured_errors_are_inside_the_derived_tolerances(
    committed: golden.SuiteResult,
) -> None:
    """Every asserted metric is inside its two-grid tolerance — and the numbers are real."""
    worst = max(
        (abs(m.relative_error) / m.tolerance)
        for case in committed.cases
        for m in case.metrics
        if m.kind == "analytic" and m.relative_error is not None and m.tolerance
    )
    assert worst <= 1.0
    assert worst > 0.0  # not a suite of exact zeros only


def test_baseline_provenance_and_policy_are_complete() -> None:
    data = json.loads(BASELINE.read_text(encoding="utf-8"))
    assert data["schema"] == golden.SCHEMA
    provenance = data["provenance"]
    for key in ("generated_by", "git_commit", "python", "platform", "geometry_fingerprint"):
        assert provenance[key]
    assert "two-grid" in data["policy"]["acceptance_tolerance"]
    assert data["policy"]["drift_tolerance"]
    for name in ("wigley", "series60", "sphere_ellipsoid", "dtmb5415"):
        case = data["cases"][name]
        assert case["reference"]
        assert set(case["measurements"]) == {"coarse", "fine"}
        assert case["tolerances"]
        assert case["reference_source"]


def test_baseline_reference_values_reproduce_from_their_source() -> None:
    """No hand-copied criteria: every blessed reference is recomputable right now."""
    data = json.loads(BASELINE.read_text(encoding="utf-8"))
    for case in golden.golden_cases():
        stored = data["cases"][case.name]["reference"]
        recomputed = case.reference()
        assert set(stored) == set(recomputed)
        for key, value in recomputed.items():
            assert golden._relative(float(stored[key]) - value, value) <= 1e-9


def test_known_limitations_are_recorded_not_hidden(committed: golden.SuiteResult) -> None:
    """The 5415 open-shell displacement is disclosed with numbers, not dressed up."""
    dtmb = next(case for case in committed.cases if case.name == "dtmb5415")
    limitation = next(
        entry
        for entry in dtmb.known_limitations
        if entry["id"] == "displacement-withheld-open-shell"
    )
    assert limitation["asserted"] is False
    assert limitation["boundary_edge_count"] > 0.0
    assert limitation["displacement_volume_m3"] > 0.0
    # The number that used to be published as a volume: it is the signed tetrahedron sum
    # of an open surface, and it is negative, which no displacement can be.
    assert limitation["signed_tetra_sum_open_shell_m3"] < 0.0


# --------------------------------------------------------------------------------------
# engine: blessing and drift
# --------------------------------------------------------------------------------------
def test_bless_then_check_round_trip(tmp_path: Path) -> None:
    path = _bless(tmp_path, ProbeCase())
    data = json.loads(path.read_text(encoding="utf-8"))
    assert data["cases"]["probe"]["measurements"]["fine"]["value"]
    result = golden.run_suite(baseline_path=path, cases=[ProbeCase()])
    assert result.passed, [f"{f.case}/{f.code}: {f.detail}" for f in result.errors]


def test_drift_beyond_tolerance_is_an_error(tmp_path: Path) -> None:
    path = _bless(tmp_path, ProbeCase())
    drifted = golden.run_suite(baseline_path=path, cases=[ProbeCase(scale=1.01)])
    codes = _codes(drifted)
    assert not drifted.passed
    assert "drift" in codes, codes
    detail = next(f for f in drifted.errors if f.code == "drift").detail
    assert "drifted" in detail


def test_drift_inside_tolerance_is_not_an_error(tmp_path: Path) -> None:
    path = _bless(tmp_path, ProbeCase())
    unchanged = golden.run_suite(baseline_path=path, cases=[ProbeCase(scale=1.0 + 1e-12)])
    assert unchanged.passed, [f.detail for f in unchanged.errors]


def test_hand_edited_reference_is_rejected(tmp_path: Path) -> None:
    path = _bless(tmp_path, ProbeCase())
    data = json.loads(path.read_text(encoding="utf-8"))
    data["cases"]["probe"]["reference"]["value"] *= 1.03
    path.write_text(json.dumps(data), encoding="utf-8")
    result = golden.run_suite(baseline_path=path, cases=[ProbeCase()])
    codes = _codes(result)
    assert "reference_not_reproducible" in codes, codes
    assert not result.passed


def test_missing_baseline_fails_closed(tmp_path: Path) -> None:
    result = golden.run_suite(baseline_path=tmp_path / "absent.json", cases=[ProbeCase()])
    assert "baseline_missing" in _codes(result)
    assert not result.passed


def test_unreadable_baseline_fails_closed(tmp_path: Path) -> None:
    path = tmp_path / "broken.json"
    path.write_text("{not json", encoding="utf-8")
    result = golden.run_suite(baseline_path=path, cases=[ProbeCase()])
    assert "baseline_unreadable" in _codes(result)


def test_wrong_schema_fails_closed(tmp_path: Path) -> None:
    path = tmp_path / "schema.json"
    path.write_text(json.dumps({"schema": "other/1", "cases": {}}), encoding="utf-8")
    result = golden.run_suite(baseline_path=path, cases=[ProbeCase()])
    assert "baseline_schema" in _codes(result)
    assert not result.passed


def test_case_unknown_to_the_baseline_fails_closed(tmp_path: Path) -> None:
    path = _bless(tmp_path, ProbeCase())
    result = golden.run_suite(baseline_path=path, cases=[ProbeCase(), _OtherProbe()])
    assert "case_not_in_baseline" in _codes(result)
    assert not result.passed


def test_stale_case_in_the_baseline_fails_closed(tmp_path: Path) -> None:
    path = _bless(tmp_path, ProbeCase())
    result = golden.run_suite(baseline_path=path, cases=[])
    assert "baseline_case_unknown" in _codes(result)
    assert not result.passed


def test_reference_metric_missing_from_the_baseline_is_an_error(tmp_path: Path) -> None:
    path = _bless(tmp_path, ProbeCase())
    data = json.loads(path.read_text(encoding="utf-8"))
    data["cases"]["probe"]["reference"] = {}
    path.write_text(json.dumps(data), encoding="utf-8")
    result = golden.run_suite(baseline_path=path, cases=[ProbeCase()])
    assert "reference_not_in_baseline" in _codes(result)


def test_reference_mismatch_outside_tolerance_is_an_error(tmp_path: Path) -> None:
    """A 5 % shift of the measurement is far outside a tolerance derived from 0.22 %."""
    path = _bless(tmp_path, ProbeCase())
    result = golden.run_suite(baseline_path=path, cases=[ProbeCase(scale=1.05)])
    assert "reference_mismatch" in _codes(result)
    metric = next(m for m in result.cases[0].metrics if m.name == "value")
    assert not metric.passed
    assert metric.tolerance is not None and metric.tolerance < 0.05


def test_invariant_violation_is_an_error(tmp_path: Path) -> None:
    path = _bless(tmp_path, ProbeCase())
    result = golden.run_suite(baseline_path=path, cases=[ProbeCase(invariant_residual=1e-6)])
    assert "invariant_violated" in _codes(result)
    assert not result.passed


def test_tolerance_is_twice_the_two_grid_change(tmp_path: Path) -> None:
    path = _bless(tmp_path, ProbeCase())
    result = golden.run_suite(baseline_path=path, cases=[ProbeCase()])
    metric = next(m for m in result.cases[0].metrics if m.name == "value")
    coarse = result.cases[0].measurements["coarse"]["value"]
    fine = result.cases[0].measurements["fine"]["value"]
    reference = result.cases[0].reference["value"]
    expected = 2.0 * abs(fine - coarse) / reference
    assert metric.tolerance == pytest.approx(expected, rel=1e-12)
    assert f"{expected:.3e}" in metric.basis


def test_convergence_is_checked_not_assumed(tmp_path: Path) -> None:
    """If the fine grid moves away from the reference, the gate says so."""

    class NonConverging(golden.GoldenCase):
        name = "probe"
        title = "non-converging probe"
        reference_source = "tests.test_golden_regression.NonConverging.reference"

        def measure(self, level: str) -> dict[str, float]:
            return {"value": 1.0 + {"coarse": 0.01, "fine": 0.05}[level]}

        def reference(self) -> dict[str, float]:
            return {"value": 1.0}

    path = tmp_path / "baseline.json"
    golden.bless(path, cases=[NonConverging()], force=True)
    result = golden.run_suite(baseline_path=path, cases=[NonConverging()])
    assert "convergence_not_monotone" in _codes(result)
    assert any(f.severity == golden.WARNING for f in result.all_findings)


def test_bless_refuses_a_tree_that_fails_its_own_gate(tmp_path: Path) -> None:
    path = _bless(tmp_path, ProbeCase())
    with pytest.raises(RuntimeError, match="refusing to bless"):
        golden.bless(path, cases=[ProbeCase(scale=1.2)])
    forced = tmp_path / "forced.json"
    baseline, result = golden.bless(forced, cases=[ProbeCase(scale=1.2)], force=True)
    assert forced.is_file()
    assert baseline["cases"]["probe"]["measurements"]["fine"]["value"] > 0.0
    assert result.cases


def test_quick_mode_skips_the_fine_level(tmp_path: Path) -> None:
    path = _bless(tmp_path, ProbeCase())
    quick = golden.run_suite(baseline_path=path, cases=[ProbeCase()], quick=True)
    assert set(quick.cases[0].measurements) == {"coarse"}
    assert quick.quick
    assert quick.passed, [f.detail for f in quick.errors]


def test_case_crash_is_a_failed_gate_not_a_skip(tmp_path: Path) -> None:
    class Broken(golden.GoldenCase):
        name = "probe"
        title = "broken probe"
        reference_source = "tests.test_golden_regression.Broken.reference"

        def measure(self, level: str) -> dict[str, float]:
            raise ValueError("pipeline exploded")

        def reference(self) -> dict[str, float]:
            return {"value": 1.0}

    path = tmp_path / "baseline.json"
    golden.bless(path, cases=[Broken()], force=True)
    result = golden.run_suite(baseline_path=path, cases=[Broken()])
    assert "case_crashed" in _codes(result)
    assert not result.passed


class _OtherProbe(golden.GoldenCase):
    name = "other_probe"
    title = "second probe"
    reference_source = "tests.test_golden_regression._OtherProbe.reference"

    def measure(self, level: str) -> dict[str, float]:
        return {"value": 2.0}

    def reference(self) -> dict[str, float]:
        return {"value": 2.0}


# --------------------------------------------------------------------------------------
# the trigger surface: what makes the gate run
# --------------------------------------------------------------------------------------
def test_staged_path_selection() -> None:
    assert golden.staged_requires_golden(["tensorcad/wigley.py"])
    assert golden.staged_requires_golden(["tensorcad/geometry.py", "docs/readme.md"])
    assert golden.staged_requires_golden(["tensorcad/new_module.py"])
    assert golden.staged_requires_golden(["docs/golden/tensorcad_geometry_golden.json"])
    assert not golden.staged_requires_golden(["docs/baseline_fleet.md", "README.md"])
    assert not golden.staged_requires_golden([])


def test_geometry_files_cover_the_package() -> None:
    relative = {path.relative_to(ROOT).as_posix() for path in golden.geometry_files()}
    assert "tensorcad/geometry.py" in relative
    assert "tensorcad/wigley.py" in relative
    assert "tensorcad/series60.py" in relative
    assert "tensorcad/dtmb5415.py" in relative
    assert "tensorcad/golden.py" in relative
    assert not any(path.startswith("tests/") and "golden" not in path for path in relative)


def test_geometry_fingerprint_is_content_addressed(tmp_path: Path) -> None:
    fake = tmp_path / "repo"
    (fake / "tensorcad").mkdir(parents=True)
    (fake / "tests").mkdir()
    (fake / "tensorcad" / "a.py").write_text("x = 1\n", encoding="utf-8")
    first = golden.geometry_fingerprint(fake)
    (fake / "tensorcad" / "a.py").write_text("x = 2\n", encoding="utf-8")
    second = golden.geometry_fingerprint(fake)
    assert first != second
    assert len(first) == 64


def test_git_hook_skips_when_nothing_geometric_is_staged(tmp_path: Path) -> None:
    listing = tmp_path / "staged.txt"
    listing.write_text("docs/baseline_fleet.md\nREADME.md\n", encoding="utf-8")
    completed = subprocess.run(
        ["/bin/sh", str(ROOT / ".githooks" / "pre-commit")],
        cwd=ROOT,
        env={
            "PATH": "/usr/bin:/bin",
            "PYTHON": sys.executable,
            "TENSORCAD_GOLDEN_STAGED_FILE": str(listing),
        },
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    assert "skipped" in completed.stdout


def test_git_hook_runs_the_gate_on_a_geometry_change(tmp_path: Path) -> None:
    """The hook path a developer actually hits: a staged geometry file runs the gate."""
    listing = tmp_path / "staged.txt"
    listing.write_text("tensorcad/geometry.py\n", encoding="utf-8")
    completed = subprocess.run(
        ["/bin/sh", str(ROOT / ".githooks" / "pre-commit")],
        cwd=ROOT,
        env={
            "PATH": "/usr/bin:/bin",
            "PYTHON": sys.executable,
            "TENSORCAD_GOLDEN_STAGED_FILE": str(listing),
        },
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert "golden geometry regression" in completed.stdout
    assert "PASS" in completed.stdout


def test_git_hook_honours_the_skip_switch(tmp_path: Path) -> None:
    listing = tmp_path / "staged.txt"
    listing.write_text("tensorcad/geometry.py\n", encoding="utf-8")
    completed = subprocess.run(
        ["/bin/sh", str(ROOT / ".githooks" / "pre-commit")],
        cwd=ROOT,
        env={
            "PATH": "/usr/bin:/bin",
            "PYTHON": sys.executable,
            "TENSORCAD_GOLDEN_STAGED_FILE": str(listing),
            "TENSORCAD_GOLDEN_SKIP": "1",
        },
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0
    assert "skipped" in completed.stdout


def test_ci_workflow_is_wired_to_geometry_paths() -> None:
    workflow = (ROOT / ".github" / "workflows" / "golden-geometry.yml").read_text(encoding="utf-8")
    assert "tensorcad/**" in workflow
    assert "python -m tensorcad.golden --check" in workflow
    assert "tests/test_golden_regression.py" in workflow


# --------------------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------------------
def test_cli_list_and_staged_skip(capsys: pytest.CaptureFixture[str], tmp_path: Path) -> None:
    assert golden.main(["--list"]) == 0
    listing = tmp_path / "staged.txt"
    listing.write_text("README.md\n", encoding="utf-8")
    assert golden.main(["--staged", str(listing)]) == 0
    assert "skipped" in capsys.readouterr().out


def test_cli_check_writes_a_json_report(tmp_path: Path) -> None:
    report = tmp_path / "report.json"
    exit_code = golden.main(["--check", "--quick", "--json", str(report)])
    assert exit_code == 0
    data = json.loads(report.read_text(encoding="utf-8"))
    assert data["schema"] == golden.SCHEMA
    assert data["passed"] is True
    assert data["quick"] is True
    assert len(data["cases"]) == 4


def test_cli_check_fails_on_a_stale_baseline(tmp_path: Path) -> None:
    path = _bless(tmp_path, ProbeCase())
    data = json.loads(path.read_text(encoding="utf-8"))
    data["cases"]["probe"]["measurements"]["coarse"]["value"] *= 2.0
    path.write_text(json.dumps(data), encoding="utf-8")
    assert golden.main(["--check", "--baseline", str(path)]) == 1


def test_cli_bless_and_force(tmp_path: Path) -> None:
    path = tmp_path / "cli.json"
    assert golden.main(["--bless", "--baseline", str(path)]) == 0
    assert path.is_file()
    assert golden.main(["--bless", "--baseline", str(path)]) == 0


def test_format_table_reports_every_case(committed: golden.SuiteResult) -> None:
    table = golden.format_table(committed)
    assert "PASS" in table
    for name in ("wigley", "series60", "sphere_ellipsoid", "dtmb5415"):
        assert name in table
    assert "known limitation" in table
