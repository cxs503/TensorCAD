"""Tests for the command-line entry point (the agent-facing surface)."""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

from tensorcad.cli import main

REPO_ROOT = Path(__file__).resolve().parents[1]


def test_summary_reports_measured_numbers_and_artifacts(tmp_path, capsys) -> None:
    # nx=21/nz=9 keeps the sampling error under 2 %, so the only checks reported
    # are the ones this test is about (an under-sampled grid now warns, see below).
    code = main(["--length", "12", "--beam", "2.4", "--draft", "0.9", "--nx", "21", "--nz", "9",
                 "--output", str(tmp_path / "run")])
    assert code == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["parameters"] == {
        "length": 12.0, "beam": 2.4, "draft": 0.9, "nx": 21, "nz": 9, "water_density": 1025.0,
    }
    assert payload["validation"] == {"passed": True, "errors": [], "warnings": []}
    assert payload["hydrostatics"]["displacement_volume_m3"] > 0
    assert payload["mesh"]["vertex_count"] == 19 * 17 + 18  # (nx-2)(2nz-1) + 2nz
    for key in ("mesh", "offsets", "metadata"):
        assert Path(payload["artifacts"][key]).is_file()


def test_json_flag_emits_the_full_reproducibility_record(tmp_path, capsys) -> None:
    assert main(["--nx", "9", "--nz", "5", "--output", str(tmp_path), "--density", "1000", "--json"]) == 0
    metadata = json.loads(capsys.readouterr().out)
    assert metadata["units"] == "m"
    assert metadata["hydrostatics"]["water_density_kg_per_m3"] == 1000.0
    # Both numbers are rounded for the file, so compare within the rounding.
    mass = metadata["hydrostatics"]["displacement_mass_kg"]
    assert abs(mass - metadata["hydrostatics"]["displacement_volume_m3"] * 1000.0) < 1e-3
    assert "engineering_status" in metadata
    assert metadata["discretisation"]["status"] == "measured"


def test_no_convergence_skips_the_refined_regeneration(tmp_path, capsys) -> None:
    assert main(["--nx", "9", "--nz", "5", "--output", str(tmp_path), "--no-convergence"]) == 0
    capsys.readouterr()
    metadata = json.loads((tmp_path / "metadata.json").read_text())
    assert "discretisation" not in metadata


def test_no_export_writes_nothing(tmp_path, capsys) -> None:
    assert main(["--nx", "9", "--nz", "5", "--output", str(tmp_path / "unused"), "--no-export"]) == 0
    assert not (tmp_path / "unused").exists()
    payload = json.loads(capsys.readouterr().out)
    assert payload["artifacts"] == {}


def test_invalid_parameters_exit_with_code_two_and_one_clean_message(tmp_path, capsys) -> None:
    code = main(["--length", "0", "--output", str(tmp_path)])
    captured = capsys.readouterr()
    assert code == 2
    assert "invalid model" in captured.err
    assert "length must be > 0" in captured.err
    assert captured.out == ""
    assert list(tmp_path.iterdir()) == []


def test_unusable_grid_ratio_is_a_warning_not_a_failure(tmp_path, capsys) -> None:
    assert main(["--nx", "40", "--nz", "5", "--output", str(tmp_path), "--no-convergence"]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["validation"]["passed"] is True
    assert payload["validation"]["warnings"]


def test_module_entry_point_runs_from_a_clean_interpreter(tmp_path) -> None:
    result = subprocess.run(
        [sys.executable, "-m", "tensorcad.cli", "--nx", "11", "--nz", "5", "--output", str(tmp_path / "sub")],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    payload = json.loads(result.stdout)
    assert payload["validation"]["passed"] is True
    assert (tmp_path / "sub" / "hull.obj").stat().st_size > 1000


def test_a_grid_that_cannot_be_refined_still_exports_a_complete_artifact_set(
    tmp_path, capsys
) -> None:
    """The convergence estimate is reporting, not generation: it cannot fail a model.

    ``nx=2002, nz=3`` builds a legal hull but cannot be refined (4003 stations), so
    the estimate reports "skipped" and the export has all three artifacts.
    """
    target = tmp_path / "wide"
    assert main(["--nx", "2002", "--nz", "3", "--output", str(target), "--json"]) == 0
    metadata = json.loads(capsys.readouterr().out)
    assert metadata["discretisation"]["status"] == "skipped"
    assert "parameter limits" in metadata["discretisation"]["reason"]
    for name in ("hull.obj", "offsets.csv", "metadata.json"):
        assert (target / name).stat().st_size > 0
    on_disk = json.loads((target / "metadata.json").read_text())
    assert on_disk["parameters"]["nx"] == 2002


def test_a_failure_after_validation_is_structured_json_not_a_traceback(
    tmp_path, capsys, monkeypatch
) -> None:
    import tensorcad.cli as cli

    def boom(*_args, **_kwargs):
        raise RuntimeError("disk full")

    monkeypatch.setattr(cli, "export_hull", boom)
    code = main(["--nx", "11", "--nz", "5", "--output", str(tmp_path), "--json"])
    captured = capsys.readouterr()
    assert code == 2
    payload = json.loads(captured.out)
    assert payload["stage"] == "export"
    assert "disk full" in payload["error"]
    assert "export failed" in captured.err
    assert "Traceback" not in captured.err
