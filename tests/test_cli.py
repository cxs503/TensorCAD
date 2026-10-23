"""Tests for the command-line entry point (the agent-facing surface)."""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

from tensorcad.cli import main

REPO_ROOT = Path(__file__).resolve().parents[1]


def test_summary_reports_measured_numbers_and_artifacts(tmp_path, capsys) -> None:
    code = main(["--length", "12", "--beam", "2.4", "--draft", "0.9", "--nx", "11", "--nz", "5",
                 "--output", str(tmp_path / "run")])
    assert code == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["parameters"] == {
        "length": 12.0, "beam": 2.4, "draft": 0.9, "nx": 11, "nz": 5, "water_density": 1025.0,
    }
    assert payload["validation"] == {"passed": True, "errors": [], "warnings": []}
    assert payload["hydrostatics"]["displacement_volume_m3"] > 0
    assert payload["mesh"]["vertex_count"] == 9 * 9 + 10  # (nx-2)(2nz-1) + 2nz
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
