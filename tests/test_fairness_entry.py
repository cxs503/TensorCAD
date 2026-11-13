"""Entry, integration and evidence closure for GAP-CAD-A3-CURVATURE-FAIRNESS.

The capability itself (``tensorcad.fairness_energy``) is tested by
``tests/test_fairness_energy.py``, which asserts the physics.  This file asserts the
*delivery surface*, because a library nobody can reach is not a capability:

* the public package exports it (``from tensorcad import smoothing_energy``),
* the console script reaches it (``tensorcad-hull fairness ...``),
* the main hull-audit pipeline calls it (``audit_hull_geometry`` ->
  ``smoothing_energy``), so a ship model that passes through the audit carries the
  number,
* the CLI writes a machine record whose claims are the acceptance clause, and the
  benchmark driver's gates agree with the analytic references.
"""
from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import tensorcad
from tensorcad.fairness_cli import main as fairness_main
from tensorcad.fairness_energy import fairness_energy_report, smoothing_energy
from tensorcad.geometry_audit import audit_hull_geometry
from tensorcad.wigley import WigleyParameters, analytic_flank_energy, build_mesh

PUBLIC_NAMES = (
    "ENERGY_KINDS",
    "FairnessEnergyReport",
    "MeshCurvatureField",
    "SmoothingEnergy",
    "fairness_energy_report",
    "principal_curvature_field",
    "smoothing_energy",
)


def test_public_package_exports_the_capability():
    for name in PUBLIC_NAMES:
        assert name in tensorcad.__all__, f"{name} missing from tensorcad.__all__"
        assert hasattr(tensorcad, name), f"{name} not importable from the package root"
    # The re-export must be the same object, not a look-alike.
    assert tensorcad.smoothing_energy is smoothing_energy


@pytest.mark.parametrize("case", ["sphere", "wigley", "hull"])
def test_cli_case_runs_and_writes_a_record(case, tmp_path):
    out = tmp_path / f"{case}.json"
    expected_returncode = 0 if case != "wigley" else 0
    code = fairness_main([case, "--out", str(out), "--quiet"])
    assert code == expected_returncode
    record = json.loads(out.read_text())
    assert record["schema"] == "tensorcad.fairness_cli/1"
    assert record["case"] == case
    assert record["claims"], "a record with no claims cannot be judged"
    assert all(row["holds"] for row in record["claims"])
    bending = record["observed"]["energies"]["bending"]
    assert math.isfinite(bending["smooth_energy"])
    assert 0.0 <= bending["crease_area_fraction"] <= 1.0


def test_cli_sphere_claim_is_the_acceptance_clause(tmp_path):
    out = tmp_path / "sphere.json"
    assert fairness_main(["sphere", "--radius", "1", "--level", "16", "--probes", "64",
                          "--out", str(out), "--quiet"]) == 0
    record = json.loads(out.read_text())
    curvature = [row for row in record["claims"] if "1/R" in row["claim"]][0]
    assert curvature["limit"] == 1e-6          # A3 clause, absolute 1/m
    assert curvature["value"] <= 1e-6
    exact = record["analytic_reference"]
    assert exact["E_bending"] == pytest.approx(8.0 * math.pi)
    assert exact["E_willmore"] == pytest.approx(4.0 * math.pi)


def test_package_console_script_reaches_the_entry():
    """``tensorcad.cli fairness`` is the product entry, not only a module."""
    completed = subprocess.run(
        [sys.executable, "-m", "tensorcad.cli", "fairness", "wigley", "--nx", "81", "--nz", "25",
         "--quiet"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        timeout=300,
    )
    assert completed.returncode == 0, completed.stderr
    record = json.loads(completed.stdout)
    assert record["case"] == "wigley"
    assert record["all_claims_hold"] is True


def test_hull_audit_pipeline_calls_the_smoothing_energy():
    """The main geometry audit carries the number: no separate import required."""
    mesh, _ = build_mesh(WigleyParameters(nx=81, nz=25))
    report = audit_hull_geometry(mesh)
    assert "smoothing_energy" in report.stages_run
    assert report.smoothing_energy is not None
    payload = report.to_dict()["smoothing_energy"]
    reference = analytic_flank_energy(WigleyParameters())
    # Same surface, same quadrature: the audit's number is the library's number.
    direct = fairness_energy_report(mesh, crease_angle_degrees=45.0)
    assert payload["bending"]["smooth_energy"] == pytest.approx(direct.bending.smooth_energy)
    assert abs(payload["bending"]["smooth_energy"] / reference.bending - 1.0) < 0.15
    # The identity 4H^2 - 2K = k1^2 + k2^2 holds on the smooth cells to round-off.
    assert abs(payload["identity_gap"]) < 1e-9


def test_smoothing_energy_can_be_switched_off():
    mesh, _ = build_mesh(WigleyParameters(nx=41, nz=13))
    report = audit_hull_geometry(mesh, check_smoothing_energy=False)
    assert "smoothing_energy" not in report.stages_run
    assert report.smoothing_energy is None
    assert "smoothing_energy" not in report.to_dict()


def test_benchmark_gates_hold_at_the_declared_settings():
    """Run the verified benchmark's gate logic at reduced cost and check the rows."""
    sys.path.insert(0, str(REPO_ROOT / "benchmarks" / "verified" / "curvature_fairness"))
    import run_benchmark  # noqa: E402  (path-injected on purpose)

    result = run_benchmark.run(
        resolutions=((21, 7), (41, 13), (81, 25)), sphere_level=16, probes=32
    )
    # The gated blocks all pass; the diagnostics are allowed to report a weakness.
    assert result["all_gates_pass"] is True, {
        name: block["failing"] for name, block in result["gates"].items()
        if block["verdict"] != "pass"
    }
    sphere_rows = result["gates"]["G1 sphere principal curvature == 1/R (A3 clause, 1e-6 1/m)"]["rows"]
    assert len(sphere_rows) == 7
    assert max(row["error"] for row in sphere_rows) < 1e-6
    assert "wigley_flank" in result["cross_check"]["three_case_agreement"]
    assert result["cross_check"]["provenance"]["status"].startswith("meshmagick comparison WITHHELD")