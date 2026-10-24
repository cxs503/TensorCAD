"""CLI tests for ``tensorcad-neural``: JSON output, exit codes and artifacts."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

pytest.importorskip("torch", reason="the neural geometry CLI needs the 'neural' extra")

from tensorcad.cli_neural import main
from tensorcad.neural.objects import Tolerances, extract_object, fit_sphere

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


def fitted(tmp_path: Path) -> Path:
    target = tmp_path / "sphere.ngo"
    fit_sphere(target, radius=0.5, steps=20, seed=5, tolerance=LOOSE, rays=48)
    return target


def parse(stdout: str) -> dict:
    return json.loads(stdout, parse_constant=lambda token: pytest.fail(f"non-standard token {token}"))


def test_info_reports_manifest_and_history(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    target = fitted(tmp_path)
    assert main(["--json", "info", "--object", str(target)]) == 0
    payload = parse(capsys.readouterr().out)
    assert payload["manifest"]["revision"] == 1
    assert payload["history"][0]["status"] == "committed"


def test_info_fails_cleanly_for_a_missing_object(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["--json", "info", "--object", str(tmp_path / "absent.ngo")]) == 2
    payload = parse(capsys.readouterr().out)
    assert payload["status"] == "error"


def test_verify_passes_on_a_published_object(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    target = fitted(tmp_path)
    assert main(["--json", "verify", "--object", str(target)]) == 0
    payload = parse(capsys.readouterr().out)
    assert payload["status"] == "verified"
    assert all(check["passed"] for check in payload["checks"])


def test_verify_detects_tampered_weights(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    target = fitted(tmp_path)
    weights = target / "weights.safetensors"
    weights.write_bytes(weights.read_bytes() + b"tampered")
    assert main(["--json", "verify", "--object", str(target)]) == 1
    payload = parse(capsys.readouterr().out)
    assert payload["status"] == "failed"
    assert any(check["name"] == "weights_hash_matches" and not check["passed"] for check in payload["checks"])


def test_extract_writes_obj_and_stl(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    target = fitted(tmp_path)
    output = tmp_path / "derived"
    code = main(
        ["--json", "extract", "--object", str(target), "--resolution", "20", "--output", str(output)]
    )
    assert code == 0
    payload = parse(capsys.readouterr().out)
    assert payload["report"]["watertight"] is True
    assert {item["format"] for item in payload["artifacts"]} == {"obj", "stl"}
    for name in ("surface.obj", "surface.stl"):
        assert (output / name).stat().st_size > 0
    assert (output / "surface.stl").read_bytes()[80:84] != b"\x00\x00\x00\x00"


def test_rejected_edit_exits_non_zero_with_json(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    target = fitted(tmp_path)
    code = main(
        [
            "--json",
            "edit",
            "--object",
            str(target),
            "--depth",
            "0.30",
            "--sigma",
            "0.08",
            "--steps",
            "20",
            "--seed",
            "3",
            "--rays",
            "48",
        ]
    )
    assert code == 1
    payload = parse(capsys.readouterr().out)
    assert payload["status"] == "rolled_back"
    assert payload["active_revision"] == 1
    assert payload["failed_checks"]


def test_fit_refuses_to_reuse_an_existing_bundle(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    target = fitted(tmp_path)
    assert main(["--json", "fit", "--output", str(target), "--steps", "5", "--rays", "8"]) == 2
    payload = parse(capsys.readouterr().out)
    assert payload["status"] == "error"
    assert "already holds revision" in payload["error"]


def test_verify_skips_the_field_query_when_the_neural_extra_is_missing(
    tmp_path: Path, capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch
) -> None:
    """Integrity checks are standard-library work and must survive a torch-free env.

    Only the live field query needs torch; the manifest, the hash, the history and
    the derived artifacts do not, so "skipped" is the correct answer rather than a
    traceback.
    """
    import sys
    import types

    target = fitted(tmp_path)
    extract_object(target, resolution=16)

    def unavailable(*_args: object, **_kwargs: object) -> object:
        raise ImportError("No module named 'torch'")

    fake = types.ModuleType("tensorcad.neural.objects")
    fake.build_field = unavailable  # type: ignore[attr-defined]
    monkeypatch.setitem(sys.modules, "tensorcad.neural.objects", fake)

    assert main(["--json", "verify", "--object", str(target)]) == 0
    payload = parse(capsys.readouterr().out)
    assert payload["status"] == "verified"
    by_name = {check["name"]: check for check in payload["checks"]}
    assert by_name["field_query"]["severity"] == "skipped"
    assert "neural extra" in by_name["field_query"]["detail"]
    assert any(name.startswith("derived_hash:") for name in by_name)
    assert payload["warnings"] == []


def test_verify_reports_a_missing_derived_artifact_as_a_warning(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    target = fitted(tmp_path)
    extract_object(target, resolution=16)
    victim = next((target / "derived").iterdir())
    victim.unlink()

    assert main(["--json", "verify", "--object", str(target)]) == 0
    payload = parse(capsys.readouterr().out)
    assert payload["status"] == "verified"
    assert payload["warnings"], "a lost export is worth reporting even if geometry is intact"
