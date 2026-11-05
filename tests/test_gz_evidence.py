"""The GZ evidence in this repository is regenerable, so regenerate it and compare.

``examples/gz/run_example.py`` is the only producer of the numbers quoted for
``tensorcad/gz.py``; ``examples/gz/outputs/`` holds a committed copy so a reader can
check a claim without running anything. The two can drift apart silently -- a checker
whose published numbers were typed by hand can never fail -- so these tests run the
example again and hold the committed files to it, and separately re-render the
published tables from the committed JSON.

What is deliberately *not* asserted: byte identity of ``plot.png`` (a matplotlib
version may change the encoder) and byte identity of ``results.json`` (round-off in
``sin``/``cos`` is libm's business). The numeric comparison below allows round-off
noise and nothing else, so an actual change of method or parameter fails it.
"""

from __future__ import annotations

import importlib.util
import json
import math
import subprocess
import sys
from pathlib import Path
from types import ModuleType

import pytest

ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "examples" / "gz" / "run_example.py"
OUTPUTS = ROOT / "examples" / "gz" / "outputs"

def _load_example() -> ModuleType:
    spec = importlib.util.spec_from_file_location("tensorcad_example_gz", EXAMPLE)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def fresh_run(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Run the example once, as CI does, into a scratch directory."""
    output = tmp_path_factory.mktemp("gz_evidence")
    result = subprocess.run(
        [sys.executable, str(EXAMPLE), "--output", str(output)],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    return output


def _compare(fresh: object, committed: object, path: str = "") -> None:
    """Structural equality, with a tolerance that admits round-off but not drift."""
    if isinstance(fresh, dict):
        assert isinstance(committed, dict), path
        assert list(fresh) == list(committed), path
        for key in fresh:
            _compare(fresh[key], committed[key], f"{path}/{key}")
    elif isinstance(fresh, list):
        assert isinstance(committed, list), path
        assert len(fresh) == len(committed), f"{path}: {len(fresh)} vs {len(committed)}"
        for index, (a, b) in enumerate(zip(fresh, committed, strict=True)):
            _compare(a, b, f"{path}[{index}]")
    elif isinstance(fresh, float) or isinstance(committed, float):
        assert isinstance(fresh, (int, float)) and isinstance(committed, (int, float)), path
        assert math.isclose(fresh, committed, rel_tol=1e-9, abs_tol=1e-15), (
            f"{path}: {fresh!r} vs {committed!r}"
        )
    else:
        assert fresh == committed, f"{path}: {fresh!r} vs {committed!r}"


def test_the_committed_run_passed_its_own_checks() -> None:
    """A published evidence set that records a failure is not evidence."""
    data = json.loads((OUTPUTS / "results.json").read_text())
    assert data["checks"]["failures"] == []
    assert data["checks"]["passed"] is True


def test_committed_results_reproduce_from_the_current_code(fresh_run: Path) -> None:
    """Drift guard: the published numbers are what the code computes *now*."""
    fresh = json.loads((fresh_run / "results.json").read_text())
    committed = json.loads((OUTPUTS / "results.json").read_text())
    _compare(fresh, committed)


def test_committed_report_is_reproducible(fresh_run: Path) -> None:
    assert (fresh_run / "report.md").read_text() == (OUTPUTS / "report.md").read_text()


def test_committed_tables_are_rendered_from_the_committed_json() -> None:
    """No hand-typed numbers: the tables are a function of ``results.json``."""
    render = _load_example().render_tables_markdown
    committed = json.loads((OUTPUTS / "results.json").read_text())
    assert render(committed) == (OUTPUTS / "tables.md").read_text()


def test_committed_tables_match_a_fresh_render(fresh_run: Path) -> None:
    """Same generator, same input: the checked-in tables are not an older run."""
    assert (fresh_run / "tables.md").read_text() == (OUTPUTS / "tables.md").read_text()


def test_tables_quote_the_verified_errors() -> None:
    """Spot-check the claim a reader would look for: the two closed-form agreements."""
    data = json.loads((OUTPUTS / "results.json").read_text())
    text = (OUTPUTS / "tables.md").read_text()
    assert f"{data['box_wall_sided']['worst_abs_error_m']:.2e} m" in text
    assert f"{data['log_circular_section']['worst_abs_error_m']:.2e} m" in text
    assert data["box_wall_sided"]["worst_abs_error_m"] < 1e-12
    assert data["log_circular_section"]["worst_abs_error_m"] < 1e-12


def _check_png(path: Path) -> None:
    blob = path.read_bytes()
    assert blob.startswith(b"\x89PNG\r\n\x1a\n"), path
    width = int.from_bytes(blob[16:20], "big")
    height = int.from_bytes(blob[20:24], "big")
    assert width > 200 and height > 100, (path, width, height)


def test_the_committed_plot_is_a_non_empty_png() -> None:
    """Readable without matplotlib: the committed figure is checked byte-header first."""
    _check_png(OUTPUTS / "plot.png")


def test_a_fresh_run_rewrites_the_plot(fresh_run: Path) -> None:
    """Skipped where matplotlib is absent -- the example says so and skips too."""
    pytest.importorskip("matplotlib", reason="the plot is optional: the example skips it")
    _check_png(fresh_run / "plot.png")
