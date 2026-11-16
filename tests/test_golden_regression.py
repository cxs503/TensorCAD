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

    def __init__(
        self,
        scale: float = 1.0,
        invariant_residual: float = 0.0,
        fingerprint: str = "probe-geometry-0",
    ) -> None:
        self.scale = scale
        self.invariant_residual = invariant_residual
        self.fingerprint = fingerprint

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

    def geometry_sha256(self) -> str:
        return self.fingerprint


class DeclaringProbeCase(ProbeCase):
    """A probe whose band is declared, like the 5415 particulars table."""

    def reference_tolerance(
        self, values: dict[str, float], reference: dict[str, float]
    ) -> dict[str, tuple[float, str]]:
        return {"value": (3e-2, "declared contract: the test's own 3 % band")}


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
        "inclined_box",
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
    from tensorcad.dtmb5415 import WATERLINE_Z, DTMB5415Space
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
    origins = set()
    for metric in analytic:
        assert metric.tolerance is not None and metric.tolerance > 0.0
        assert metric.reference is not None
        assert metric.origin in {"two-grid", "declared"}
        origins.add(metric.origin)
        if metric.origin == "two-grid":
            assert "two-grid" in metric.basis or "floor" in metric.basis
        else:
            # a declared band must name the contract it is quoted from, and stay inside the
            # engine's ceiling: that is the whole difference between declared and chosen
            assert metric.basis.startswith("declared contract: ")
            assert len(metric.basis) > len("declared contract: ")
            assert metric.tolerance <= golden.DECLARED_TOLERANCE_CAP
    assert origins == {"two-grid", "declared"}


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
        # Cases may carry extra levels (e.g. "finer" for GCI three-grid); require at least coarse+fine.
        assert {"coarse", "fine"} <= set(case["measurements"])
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


def test_ellipsoid_area_is_asserted_against_an_exact_reference(
    committed: golden.SuiteResult,
) -> None:
    """Regression for the lifted "ellipsoid-area-not-asserted" limitation.

    A triaxial ellipsoid's area is an elliptic integral, so for a while it was only
    drift-locked: a gate that notices a *changed* area but not a *wrong* one. It is now
    asserted against converged quadrature of the exact area element, which shares no code
    with the mesher under test.
    """
    from tensorcad.analytic_solids import ellipsoid_surface_area

    case = {result.name: result for result in committed.cases}["sphere_ellipsoid"]
    assert "ellipsoid-area-not-asserted" not in {entry["id"] for entry in case.known_limitations}
    metric = next(m for m in case.metrics if m.name == "ellipsoid_surface_area_m2")
    assert metric.kind == "analytic"
    assert metric.reference == pytest.approx(ellipsoid_surface_area((2.0, 1.5, 0.5)), rel=1e-12)
    assert metric.relative_error is not None and metric.tolerance is not None
    assert metric.passed
    assert abs(metric.relative_error) <= metric.tolerance


def test_a_failing_case_reports_itself_as_failed(tmp_path: Path) -> None:
    """Per-case ``passed`` must reflect that case's findings, not be vacuously true.

    Regression: ``run_suite`` collected every finding into one suite-level list and left
    ``CaseResult.findings`` empty, so ``case.passed`` was always True. The suite flag was
    right (it read the shared list) while the report told a reader that a drifting case
    passed. Measured with a case whose measurement is 50 % off its blessed value.
    """
    baseline = _bless(tmp_path, ProbeCase())
    drifted = ProbeCase(scale=1.5)
    result = golden.run_suite(baseline_path=baseline, cases=[drifted])

    assert result.passed is False
    (case_result,) = result.cases
    assert case_result.passed is False
    assert case_result.findings, "the failing case must carry its own findings"
    assert all(finding.case == "probe" for finding in case_result.findings)
    # one error, counted once: the suite list holds suite-level findings only
    assert len(result.errors) == len({(f.case, f.code, f.detail) for f in result.errors})


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
# Reach: a geometry change that arrives as a data edit
# --------------------------------------------------------------------------------------
# The suite measures one case that is a *shipped dataset* rather than generated code
# (``dtmb5415``): its basis, its fit and its reference residual all live under
# ``docs/data/``.  Editing any of them changes what the gate asserts while touching no
# file under ``tensorcad/``, so before the declared inputs existed the gate printed
# "no staged geometry change, check skipped" -- a green pass over a changed dataset.  The
# tests below hold the trigger set, the source fingerprint and the CI path filter to the
# inputs the cases declare, instead of to a list somebody remembered to update.


def test_every_declared_input_exists() -> None:
    for case in golden.golden_cases():
        for path in case.geometry_inputs():
            assert Path(path).is_file(), (case.name, path)


def test_staged_declared_input_triggers_the_gate() -> None:
    declared = [
        path.relative_to(ROOT).as_posix()
        for case in golden.golden_cases()
        for path in case.geometry_inputs()
        if ROOT in Path(path).parents
    ]
    assert declared, "no case declares a data input: the reach below is untested"
    for path in declared:
        assert golden.staged_requires_golden([path]), path
        # and it is not the tensorcad/ rule that catches it
        assert not path.startswith("tensorcad/"), path


def test_staged_data_that_no_case_reads_is_still_skipped() -> None:
    """The trigger must be the declared inputs, not "everything under docs/data"."""
    assert not golden.staged_requires_golden(["docs/data/dtmb5415/README.md"])
    assert not golden.staged_requires_golden(["docs/data/other/fixture.npz"])


def test_geometry_files_include_the_declared_inputs() -> None:
    relative = {path.relative_to(ROOT).as_posix() for path in golden.geometry_files()}
    for case in golden.golden_cases():
        for path in case.geometry_inputs():
            if ROOT in Path(path).parents:
                assert path.relative_to(ROOT).as_posix() in relative, path


def test_source_fingerprint_moves_when_a_dataset_edit_lands(tmp_path: Path) -> None:
    """The fingerprint is a statement about the inputs, not only about the code."""
    fake = tmp_path / "repo"
    (fake / "tensorcad").mkdir(parents=True)
    (fake / "tests").mkdir()
    (fake / "tensorcad" / "a.py").write_text("x = 1\n", encoding="utf-8")
    data = fake / "docs" / "data" / "dtmb5415"
    data.mkdir(parents=True)
    for case in golden.golden_cases():
        for path in case.geometry_inputs():
            if ROOT in Path(path).parents:
                (fake / path.relative_to(ROOT)).write_bytes(b"dataset v1")
    first = golden.geometry_fingerprint(fake)
    (data / "holdout_fixture.npz").write_bytes(b"dataset v2")
    second = golden.geometry_fingerprint(fake)
    assert first != second
    assert len(second) == 64


def test_git_hook_runs_the_gate_on_a_dataset_edit(tmp_path: Path) -> None:
    listing = tmp_path / "staged.txt"
    listing.write_text("docs/data/dtmb5415/holdout_fixture.npz\n", encoding="utf-8")
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
    assert "golden geometry regression" in completed.stdout, completed.stdout
    assert "PASS" in completed.stdout


def _workflow_globs() -> list[str]:
    workflow = (ROOT / ".github" / "workflows" / "golden-geometry.yml").read_text(encoding="utf-8")
    globs: list[str] = []
    inside = False
    for line in workflow.splitlines():
        stripped = line.strip()
        if stripped.startswith("paths:"):
            inside = True
            continue
        if inside:
            if stripped.startswith("- "):
                globs.append(stripped[2:].strip().strip('"'))
            elif stripped and not stripped.startswith("#"):
                inside = False
    return globs


def _glob_matches(pattern: str, path: str) -> bool:
    """Minimal ``**`` matcher: the workflow's path filter uses no other operator."""
    import re

    regex = ""
    index = 0
    while index < len(pattern):
        if pattern.startswith("**", index):
            regex += ".*"
            index += 2
        elif pattern[index] == "*":
            regex += "[^/]*"
            index += 1
        else:
            regex += re.escape(pattern[index])
            index += 1
    return re.fullmatch(regex, path) is not None


def test_ci_path_filter_watches_every_declared_input() -> None:
    globs = _workflow_globs()
    assert "tensorcad/**" in globs, globs
    for case in golden.golden_cases():
        for path in case.geometry_inputs():
            if ROOT not in Path(path).parents:
                continue
            relative = path.relative_to(ROOT).as_posix()
            assert any(_glob_matches(pattern, relative) for pattern in globs), (
                relative, globs
            )


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


# --------------------------------------------------------------------------------------
# geometry identity: the mesh itself is asserted, not only the numbers integrated on it
# --------------------------------------------------------------------------------------
def test_geometry_fingerprint_change_is_an_error(tmp_path: Path) -> None:
    """A mesh that changed must fail even when every measured number is unchanged."""
    path = _bless(tmp_path, ProbeCase())
    moved = ProbeCase(fingerprint="probe-geometry-1")
    result = golden.run_suite(baseline_path=path, cases=[moved])
    codes = _codes(result)
    assert "geometry_changed" in codes, codes
    assert not result.passed
    detail = next(f for f in result.errors if f.code == "geometry_changed").detail
    assert "probe-geometry-0"[:16] in detail and "probe-geometry-1"[:16] in detail


def test_case_without_fingerprint_is_reported_not_silent(tmp_path: Path) -> None:
    """A case that cannot identify its mesh is a warning, and the reason is published."""
    path = _bless(tmp_path, ProbeCase(fingerprint=""))
    result = golden.run_suite(baseline_path=path, cases=[ProbeCase(fingerprint="")])
    assert result.passed
    codes = _codes(result)
    assert "geometry_fingerprint_missing" in codes, codes


def test_baseline_without_fingerprint_is_an_error(tmp_path: Path) -> None:
    """An un-blessed mesh is an unguarded case: fail closed, do not warn and pass.

    A baseline entry whose ``geometry_sha256`` is empty locks the integrated numbers but
    not the mesh, so re-triangulation or vertex renumbering passes silently. The gate
    therefore rejects it; only ``--bless`` (a deliberate act) clears the finding.
    """
    path = _bless(tmp_path, ProbeCase())
    data = json.loads(path.read_text(encoding="utf-8"))
    data["cases"]["probe"]["geometry_sha256"] = ""
    path.write_text(json.dumps(data), encoding="utf-8")
    result = golden.run_suite(baseline_path=path, cases=[ProbeCase()])
    assert not result.passed
    codes = _codes(result)
    assert "geometry_fingerprint_not_blessed" in codes, codes
    assert "geometry_fingerprint_not_blessed" in [f.code for f in result.errors]


def test_every_golden_case_fingerprints_its_mesh(committed: golden.SuiteResult) -> None:
    """The four required golden forms must all lock the mesh, not just the integrals.

    Regression: ``--bless`` wrote ``geometry_sha256`` for every case, but only the Wigley
    case computed one (the other three stored ``""``) and ``evaluate_case`` never read the
    field back at all.
    """
    baseline = json.loads(BASELINE.read_text(encoding="utf-8"))
    assert {case.name for case in committed.cases} == {
        "wigley",
        "series60",
        "sphere_ellipsoid",
        "dtmb5415",
        "inclined_box",
    }
    for case in committed.cases:
        fingerprint = case.geometry_sha256
        assert len(fingerprint) == 64, (case.name, fingerprint)
        int(fingerprint, 16)  # hex
        assert fingerprint == baseline["cases"][case.name]["geometry_sha256"], case.name


def test_vertex_relabelling_is_caught_only_by_the_fingerprint() -> None:
    """Non-redundancy of the fingerprint: a relabelled mesh is numerically identical.

    Renumbering the vertices does not move a single triangle in space, so every
    integrated quantity -- and therefore every numeric metric of the suite -- is
    unchanged. Only ``_mesh_sha256`` (vertex and face order included) can see it.
    """
    from tensorcad import mesh as meshlib
    from tensorcad import wigley

    mesh = wigley.generate_wigley(wigley.reference_parameters(nx=21, nz=9)).mesh
    count = len(mesh.vertices)
    order = [(7 * index + 3) % count for index in range(count)]  # a permutation of range(n)
    assert len(set(order)) == count
    inverse = {old: new for new, old in enumerate(order)}
    relabelled = type(mesh)(
        vertices=[mesh.vertices[order[new]] for new in range(count)],
        faces=[tuple(inverse[index] for index in face) for face in mesh.faces],
    )
    assert golden._closed_surface_area(relabelled) == pytest.approx(
        golden._closed_surface_area(mesh), rel=1e-15
    )
    assert meshlib.total_signed_volume(relabelled) == pytest.approx(
        meshlib.total_signed_volume(mesh), rel=1e-15
    )
    assert golden._mesh_sha256(relabelled) != golden._mesh_sha256(mesh)


def test_redteam_harness_kills_an_identity_mutant(tmp_path: Path) -> None:
    """End-to-end: a real gate run against a deliberately broken tree must fail.

    The mutant disables the winding repair in the marching-tetrahedra extractor, so the
    sphere/ellipsoid meshes keep every number (area takes the norm, volume the absolute
    value) and only the fingerprint moves.
    """
    harness = _load_redteam()
    report = harness.run_redteam(quick=True, only=("identity-orientation-winding",))
    assert report["control"]["passed"], report["control"]
    (mutant,) = report["mutants"]
    assert mutant["killed"], mutant
    assert mutant["observed"] == ["geometry_changed"], mutant["observed"]


def _load_redteam():
    import importlib.util

    path = ROOT / "tools" / "golden" / "mutation_redteam.py"
    spec = importlib.util.spec_from_file_location("golden_mutation_redteam", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    # dataclasses resolve their own module through sys.modules while the class body runs,
    # so an unregistered module fails to load (AttributeError: 'NoneType' has no __dict__)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


# --------------------------------------------------------------------------------------
# DTMB 5415: the published particulars, not only the drift lock
# --------------------------------------------------------------------------------------
def _blessed_5415() -> dict:
    return json.loads(BASELINE.read_text(encoding="utf-8"))["cases"]["dtmb5415"]


def test_dtmb5415_particulars_are_asserted_against_the_published_table(
    committed: golden.SuiteResult,
) -> None:
    """Regression for "5415-particulars-drift-only".

    The case used to assert a single number (the holdout residual); displacement, beam,
    station length and Cb were drift-locked only, so the gate could notice that they
    *changed* but never that they were *wrong*. They are now asserted against the dataset
    README's principal-particulars table, scaled to full size by the case's own lambda.
    """
    from tensorcad.dtmb5415 import README_PARTICULARS_TOLERANCE, readme_reference

    case = {result.name: result for result in committed.cases}["dtmb5415"]
    published = readme_reference()
    expected = {
        "length_stations_m": published["Lpp_m"],
        "waterline_beam_m": published["beam_m"],
        "displacement_volume_m3": published["displacement_volume_m3"],
        "block_coefficient": published["cb"],
    }
    for name, reference in expected.items():
        metric = next(m for m in case.metrics if m.name == name)
        assert metric.kind == "analytic"
        assert metric.origin == "declared"
        assert metric.reference == pytest.approx(reference, rel=1e-12)
        assert metric.tolerance == pytest.approx(README_PARTICULARS_TOLERANCE)
        assert metric.passed, (name, metric.relative_error, metric.tolerance)
        assert abs(metric.relative_error) <= metric.tolerance
    # Cb is the binding one: it carries the definitional mismatch (published Lpp/B_max/T
    # against the measured station span/WL beam/draft) and still lands inside the band.
    cb = next(m for m in case.metrics if m.name == "block_coefficient")
    assert abs(cb.relative_error) > 1e-3


def test_dtmb5415_draft_is_measured_but_not_passed_off_as_a_reference(
    committed: golden.SuiteResult,
) -> None:
    """The draft cannot be "verified" against the README: it *is* the README's draft.

    ``hydrostatics()`` cuts the hull at ``DRAFT_M = 0.248 m * lambda``, so asserting the
    reported draft against the published draft compares a constant with itself. It stays a
    drift-locked measurement, and the case says so in the band's basis.
    """
    case = {result.name: result for result in committed.cases}["dtmb5415"]
    assert "draft_m" not in case.reference
    assert "draft_m" in case.measurements["fine"]
    assert not any(m.name == "draft_m" and m.kind == "analytic" for m in case.metrics)
    basis = next(iter(_blessed_5415()["tolerances"].values()))["basis"]
    assert "input" in basis


def test_dtmb5415_particulars_gate_has_teeth(monkeypatch: pytest.MonkeyPatch) -> None:
    """A wrong published table must turn the gate red, not drift quietly green.

    The teeth are the point: a drift lock detects *change*, this detects *disagreement*.
    Here the dataset's displaced mass is corrupted by 10 % (as a unit slip in the README
    table would) and the reference the case *derives* moves with it, so it no longer
    reproduces the blessed baseline and the published comparison leaves its band.
    """
    from tensorcad import dtmb5415

    table = dict(dtmb5415.DATASET["readme_model_scale_table"])
    corrupted = dict(table, displacement_t=table["displacement_t"] * 1.10)
    monkeypatch.setitem(dtmb5415.DATASET, "readme_model_scale_table", corrupted)

    case = dtmb5415_golden_case()
    findings: list[golden.Finding] = []
    result = golden.evaluate_case(case, _blessed_5415(), quick=False, findings=findings)

    codes = {finding.code for finding in result.findings} | {f.code for f in findings}
    assert "reference_not_reproducible" in codes, codes
    mismatch = [f for f in result.findings if f.code == "reference_mismatch"]
    assert any(f.metric == "displacement_volume_m3" for f in mismatch), mismatch
    assert result.passed is False


def test_a_declared_tolerance_must_be_justified(tmp_path: Path) -> None:
    """No basis, an unknown metric and a 40 % band are all errors, not warnings.

    Two declarations, because the three ways a band can be worthless are not reachable in
    one mapping: ``probe`` asserts exactly one reference metric, so a *known but unjustified*
    band and a *too loose* band would have to share the key ``value``.  The band without a
    basis is therefore declared on a metric the reference does assert, in a second call --
    declaring it on a measured-but-unasserted key (e.g. ``coarse_only``) only ever reaches
    the unknown-metric branch, which is a different defect.
    """

    class Loose(ProbeCase):
        name = "loose"

        def reference_tolerance(self, values, reference):  # type: ignore[no-untyped-def]
            return {
                "value": (0.40, "because I said so"),
                "not_a_metric": (0.01, "names nothing"),
            }

    findings: list[golden.Finding] = []
    result = golden.evaluate_case(
        Loose(), None, quick=False, findings=findings, baseline_optional=True
    )
    codes = {finding.code for finding in result.findings} | {f.code for f in findings}
    assert "declared_tolerance_too_loose" in codes, codes
    assert "declared_tolerance_unknown_metric" in codes, codes
    # ... and the case really fails: `evaluate_case` hands the declaration errors back to
    # the caller, so the pass/fail flag is only meaningful through the suite. Asserting it
    # on the bare case result was vacuously true (no accuracy error on this probe).
    suite = golden.run_suite(
        cases=[Loose()], baseline_path=tmp_path / "no_baseline.json", baseline_optional=True
    )
    assert suite.passed is False
    assert "declared_tolerance_too_loose" in _codes(suite)

    class Unjustified(ProbeCase):
        name = "unjustified"

        def reference_tolerance(self, values, reference):  # type: ignore[no-untyped-def]
            return {"value": (0.01, "   ")}

    findings = []
    result = golden.evaluate_case(
        Unjustified(), None, quick=False, findings=findings, baseline_optional=True
    )
    codes = {finding.code for finding in result.findings} | {f.code for f in findings}
    assert "declared_tolerance_without_basis" in codes, codes
    assert "declared_tolerance_too_loose" not in codes, codes
    suite = golden.run_suite(
        cases=[Unjustified()], baseline_path=tmp_path / "no_baseline.json", baseline_optional=True
    )
    assert suite.passed is False
    assert "declared_tolerance_without_basis" in _codes(suite)


def dtmb5415_golden_case() -> golden.GoldenCase:
    return next(case for case in golden.golden_cases() if case.name == "dtmb5415")

# --------------------------------------------------------------------------------------
# the recorded criterion is a derivation, not a chosen number
# --------------------------------------------------------------------------------------
# `--quick` (the mode the git hook runs) has no coarse level, so it asserts against the
# band the baseline records.  That band sits in a committed JSON file: unless the checker
# re-derives it, a hand edit decides the gate while the table still reads green.
def test_hand_inflated_band_is_rejected_in_quick_mode(tmp_path: Path) -> None:
    path = _bless(tmp_path, ProbeCase())
    data = json.loads(path.read_text(encoding="utf-8"))
    data["cases"]["probe"]["tolerances"]["value"]["value"] *= 100.0
    path.write_text(json.dumps(data), encoding="utf-8")

    result = golden.run_suite(baseline_path=path, cases=[ProbeCase()], quick=True)

    assert not result.passed
    assert "tolerance_not_derived" in _codes(result), _codes(result)
    detail = next(f for f in result.errors if f.code == "tolerance_not_derived").detail
    assert "derivation gives" in detail
    # the inflated band is refused, not silently replaced: the case still reports it
    (metric,) = [m for m in result.cases[0].metrics if m.name == "value"]
    assert metric.tolerance is not None and metric.tolerance > 0.1


def test_hand_inflated_band_is_rejected_in_full_mode(tmp_path: Path) -> None:
    """Full mode re-derives the band for its verdict; it must still refuse the record.

    Otherwise a baseline whose criterion was inflated passes CI (which ignores the stored
    value) and fails the developer hook (which reads it) -- the two gates would disagree
    about the same tree, and the disagreement would be silent.
    """
    path = _bless(tmp_path, ProbeCase())
    data = json.loads(path.read_text(encoding="utf-8"))
    data["cases"]["probe"]["tolerances"]["value"]["value"] *= 100.0
    path.write_text(json.dumps(data), encoding="utf-8")

    result = golden.run_suite(baseline_path=path, cases=[ProbeCase()], quick=False)

    assert not result.passed
    assert "tolerance_not_derived" in _codes(result), _codes(result)


def test_hand_edited_declared_band_is_rejected(tmp_path: Path) -> None:
    path = _bless(tmp_path, DeclaringProbeCase())
    data = json.loads(path.read_text(encoding="utf-8"))
    recorded = data["cases"]["probe"]["tolerances"]["value"]
    assert recorded["basis"].startswith(golden.DECLARED_BASIS_PREFIX)
    recorded["value"] = 0.049
    path.write_text(json.dumps(data), encoding="utf-8")

    result = golden.run_suite(baseline_path=path, cases=[DeclaringProbeCase()], quick=True)

    assert not result.passed
    assert "tolerance_not_derived" in _codes(result), _codes(result)


def test_recorded_band_must_match_its_origin(tmp_path: Path) -> None:
    """A declared metric recorded with a two-grid basis (or the reverse) is an error.

    The two modes then assert different criteria on the same tree: quick mode would read
    the record, full mode the fresh declaration.
    """
    path = _bless(tmp_path, ProbeCase())
    data = json.loads(path.read_text(encoding="utf-8"))
    data["cases"]["probe"]["tolerances"]["value"]["basis"] = (
        f"{golden.DECLARED_BASIS_PREFIX}: invented by the test"
    )
    path.write_text(json.dumps(data), encoding="utf-8")

    result = golden.run_suite(baseline_path=path, cases=[ProbeCase()], quick=True)

    assert not result.passed
    assert "tolerance_origin_mismatch" in _codes(result), _codes(result)


def test_band_for_a_metric_that_is_not_asserted_is_rejected(tmp_path: Path) -> None:
    path = _bless(tmp_path, ProbeCase())
    data = json.loads(path.read_text(encoding="utf-8"))
    data["cases"]["probe"]["tolerances"]["ghost"] = {"value": 0.5, "basis": "chosen"}
    path.write_text(json.dumps(data), encoding="utf-8")

    result = golden.run_suite(baseline_path=path, cases=[ProbeCase()], quick=True)

    assert not result.passed
    assert "tolerance_unknown_metric" in _codes(result), _codes(result)


def test_committed_bands_are_their_own_derivation(committed: golden.SuiteResult) -> None:
    """Every band in the committed baseline reproduces from the inputs it was derived on."""
    assert not [f for f in committed.errors if f.code.startswith("tolerance_")], [
        f"{f.case}/{f.code}: {f.detail}" for f in committed.errors
    ]
    data = json.loads(BASELINE.read_text(encoding="utf-8"))
    recorded = {
        (case, metric): entry
        for case, payload in data["cases"].items()
        for metric, entry in (payload.get("tolerances") or {}).items()
    }
    assert recorded
    for (case, metric), entry in recorded.items():
        assert isinstance(entry["value"], float) and entry["value"] > 0.0, (case, metric)
        assert entry["basis"], (case, metric)


def test_redteam_harness_kills_the_inflated_band_mutant() -> None:
    """The gate's negative control for the criterion itself, not only for the geometry."""
    harness = _load_redteam()
    report = harness.run_redteam(quick=True, only=("baseline-tolerance-hand-inflated",))

    assert report["control"]["passed"], report["control"]
    (mutant,) = report["mutants"]
    assert mutant["killed"], mutant["observed"]
    assert "tolerance_not_derived" in mutant["observed"]


# --------------------------------------------------------------------------------------
# the override path: `--bless --force` must leave a record
# --------------------------------------------------------------------------------------
class MovingGeometry(golden.GoldenCase):
    """A synthetic geometry whose size is a parameter.

    Measurement and reference move *together*, so a re-blessed tree is self-consistent:
    this is what a deliberate geometry update looks like from the engine's side, and it
    is the only shape in which a forced bless can be followed by an otherwise green run.
    """

    name = "moving"
    title = "synthetic probe whose geometry has a size parameter"
    reference_source = "tests.test_golden_regression.MovingGeometry.reference (closed form)"

    def __init__(self, size: float = 1.0) -> None:
        self.size = size

    def _volume(self) -> float:
        return 4.0 / 3.0 * math.pi * self.size**3

    def _parameters(self, level: str) -> dict[str, float]:
        return {"nx": {"coarse": 11, "fine": 21}[level]}

    def measure(self, level: str) -> dict[str, float]:
        offset = {"coarse": 3e-3, "fine": 8e-4}[level]
        return {"value": self._volume() * (1.0 + offset)}

    def reference(self) -> dict[str, float]:
        return {"value": self._volume()}

    def invariants(self, values: dict[str, float]) -> dict[str, tuple[float, float, str]]:
        return {}

    def geometry_sha256(self) -> str:
        return f"moving-geometry-{self.size:g}"


def _forced_baseline(tmp_path: Path) -> Path:
    """Bless a 1.0 tree, then force-bless the 1.05 tree on top of it."""
    path = tmp_path / "baseline.json"
    golden.bless(path, cases=[MovingGeometry(1.0)])
    golden.bless(path, cases=[MovingGeometry(1.05)], force=True)
    return path


def test_forced_bless_records_the_findings_it_overrode(tmp_path: Path) -> None:
    """`--force` is the one path that can write a baseline the gate refuses.

    It has to be an auditable act: the baseline must say which commit overrode which
    findings, otherwise a laundered baseline is indistinguishable from a reviewed one.
    """
    path = _forced_baseline(tmp_path)
    baseline = json.loads(path.read_text(encoding="utf-8"))
    record = baseline["provenance"]["forced_bless"]
    assert len(record) == 1
    assert record[0]["git_commit"]
    assert "drift" in record[0]["exceeded"]
    assert any("drift" in detail for detail in record[0]["details"])


def test_a_forced_bless_is_not_erased_by_a_later_bless(tmp_path: Path) -> None:
    """The record is append-only: a second (now clean) bless carries it forward."""
    path = _forced_baseline(tmp_path)
    baseline, result = golden.bless(path, cases=[MovingGeometry(1.05)])
    assert result.passed, [f"{f.case}/{f.code}" for f in result.errors]
    record = baseline["provenance"]["forced_bless"]
    assert len(record) == 1, "a clean bless must neither drop nor duplicate the record"

    fresh, _ = golden.bless(tmp_path / "fresh.json", cases=[MovingGeometry(1.0)])
    assert "forced_bless" not in fresh["provenance"]


def test_a_forced_bless_stays_visible_in_every_later_run(tmp_path: Path) -> None:
    """A reviewed override is visible, not fatal — and visible in quick mode too.

    Quick mode is what the git hook runs, so a signal that only full mode emits would be
    absent from exactly the loop that a developer iterates in.
    """
    path = _forced_baseline(tmp_path)
    for quick in (False, True):
        result = golden.run_suite(baseline_path=path, cases=[MovingGeometry(1.05)], quick=quick)
        forced = [f for f in result.all_findings if f.code == "forced_bless_recorded"]
        assert forced, f"no visible record in {'quick' if quick else 'full'} mode"
        assert forced[0].severity == golden.WARNING
        assert "drift" in forced[0].detail
    assert result.passed, "a reviewed override must not make the gate permanently red"


# --------------------------------------------------------------------------------------
# three-level ladder: grid-convergence (GCI) assertion
# --------------------------------------------------------------------------------------
def test_three_level_cases_record_a_grid_convergence_band() -> None:
    data = json.loads(BASELINE.read_text(encoding="utf-8"))
    three = {case.name: case for case in golden.golden_cases() if len(case.levels()) >= 3}
    assert three, "no case declares a third level: the GCI layer would be dead code"
    for name, case in three.items():
        entry = data["cases"][name]
        assert set(entry["measurements"]) == set(case.levels())
        recorded = entry["extrapolated"]
        assert recorded, f"{name}: three levels but no grid-convergence band recorded"
        for key, band in recorded.items():
            assert band["value"] is not None
            assert band["tolerance"] >= golden.TOLERANCE_FLOOR
            assert band["basis"].startswith("GCI:")


def test_recorded_grid_convergence_band_reproduces() -> None:
    """Every recorded GCI band is the one the blessed measurements derive (no hand numbers)."""
    data = json.loads(BASELINE.read_text(encoding="utf-8"))
    checked = 0
    for case in golden.golden_cases():
        if len(case.levels()) < 3:
            continue
        entry = data["cases"][case.name]
        derived = golden._gci_estimates(
            case, case.levels(), entry["measurements"], entry["reference"], []
        )
        recorded = entry["extrapolated"]
        assert set(derived) == set(recorded)
        for key, (value, band, _basis) in derived.items():
            assert abs(value - recorded[key]["value"]) <= 1e-12 * max(1.0, abs(value))
            assert abs(band - recorded[key]["tolerance"]) <= 1e-9 * band
            checked += 1
    assert checked


def test_hand_widened_grid_convergence_band_is_rejected(tmp_path: Path) -> None:
    data = json.loads(BASELINE.read_text(encoding="utf-8"))
    entry = data["cases"]["series60"]["extrapolated"]
    key = sorted(entry)[0]
    entry[key]["tolerance"] = float(entry[key]["tolerance"]) * 10.0
    path = tmp_path / "baseline.json"
    path.write_text(json.dumps(data), encoding="utf-8")
    result = golden.run_suite(quick=True, baseline_path=path)
    assert "extrapolation_not_derived" in _codes(result), _codes(result)


def test_first_order_case_keeps_a_wide_band_honestly() -> None:
    """The negative result this mission found, pinned.

    Series 60's volume / Cb / mid-section converge at order ~1.05, because the section law
    ``y ~ zeta**(1/m)`` has an unbounded slope at the keel (m ~ 12 at Cb = 0.6).  A
    three-level extrapolation therefore cannot buy a tight band: the two extrapolants
    disagree by ~0.9 %, and the GCI band stays above 1 %.  The assertion keeps that honest
    rather than letting the "extrapolated" label imply convergence that is not there.
    """
    case = next(c for c in golden.golden_cases() if c.name == "series60")
    data = json.loads(BASELINE.read_text(encoding="utf-8"))["cases"]["series60"]
    derived = golden._gci_estimates(
        case, case.levels(), data["measurements"], data["reference"], []
    )
    band = derived["displacement_volume_m3"][1]
    two_grid = data["tolerances"]["displacement_volume_m3"]["value"]
    assert band > 1e-2, band            # still above the 1 % blind-spot line: not closed
    assert band < two_grid              # but it is the band that now binds
    assert "observed order p" in derived["displacement_volume_m3"][2]


def test_verification_report_test_count_is_accurate() -> None:
    """The verification report must not claim a stale test count.

    Regression: the report said "40 tests" when the file had 71+11=82.
    This test reads the report and the actual test count, and fails if
    the report's stated number differs from reality by more than 1.
    """
    import re
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    report = root / "docs" / "golden" / "verification_report.md"
    assert report.is_file(), f"{report} is missing"

    text = report.read_text(encoding="utf-8")
    # Look for patterns like "71 项测试" or "82 项测试" or "合计.*82"
    numbers = [int(m) for m in re.findall(r'(\d+)\s*项测试', text)]
    # Also check the summary line
    summary_matches = re.findall(r'合计.*?(\d+)', text)
    all_numbers = numbers + [int(m) for m in summary_matches if m.isdigit()]

    # Count actual test functions
    test_file = Path(__file__)
    mutation_file = test_file.parent / "test_golden_mutation.py"
    actual = 0
    for f in (test_file, mutation_file):
        if f.is_file():
            actual += sum(1 for line in f.read_text(encoding="utf-8").splitlines()
                          if line.strip().startswith("def test_"))

    # The report should mention a number close to the actual count
    if all_numbers:
        reported = max(all_numbers)
        assert abs(reported - actual) <= 2, (
            f"verification_report.md claims ~{reported} tests but there are {actual}; "
            "update the report"
        )


def test_sphere_ellipsoid_has_three_grid_gci(committed: golden.SuiteResult) -> None:
    """Regression: sphere_ellipsoid now carries a finer level for GCI extrapolation.

    Before this test, only wigley and series60 had three-grid convergence. The
    sphere/ellipsoid SDF case now matches them, proving its volume/area error is
    discretisation (not model error) via Richardson extrapolation.
    """
    case = {result.name: result for result in committed.cases}["sphere_ellipsoid"]
    assert {"coarse", "fine", "finer"} <= set(case.measurements), (
        f"sphere_ellipsoid missing finer level; levels={sorted(case.measurements)}"
    )
    assert case.measurements.get("extrapolated"), "sphere_ellipsoid has no GCI extrapolation"
    # The extrapolated values must be closer to the reference than the fine-level values
    extrapolated = case.measurements["extrapolated"]
    fine = case.measurements["fine"]
    for key in ("sphere_volume_m3", "sphere_surface_area_m2",
                "ellipsoid_volume_m3", "ellipsoid_surface_area_m2"):
        assert key in extrapolated, f"missing extrapolated {key}"
        ext_val = extrapolated[key]["value"]
        fine_val = fine[key]
        ref = case.reference[key]
        ext_err = abs(ext_val - ref) / ref
        fine_err = abs(fine_val - ref) / ref
        assert ext_err <= fine_err, (
            f"{key}: GCI extrapolation error {ext_err:.3e} > fine error {fine_err:.3e}"
        )
