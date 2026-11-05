"""The published GZ evidence bundle is generated, not written: rerun it and compare.

``docs/gz/{results.json,report.md,gz_curves.png}`` is the page a reviewer reads, so it has to
be the *output* of ``examples/gz/run_example.py`` on the current tree -- otherwise the page
and the code drift apart and the numbers stop meaning anything. These tests rerun the script
into a temporary directory and hold the bundle to it:

* the JSON must reproduce, key for key, to ``1e-9`` relative;
* the Markdown must be the same page (every line identical once digits are stripped, so a
  last-digit difference between machines is tolerated but a changed claim is not);
* the figure must be a real PNG of a plot-sized canvas;
* the bundle's own version stamp must be the version of the tree that published it, which is
  what catches "bumped the version, forgot to regenerate".
"""

from __future__ import annotations

import json
import math
import re
import struct
import subprocess
import sys
from pathlib import Path

import pytest

import tensorcad

ROOT = Path(__file__).resolve().parents[1]
BUNDLE = ROOT / "docs" / "gz"
EXAMPLE = ROOT / "examples" / "gz" / "run_example.py"
DIGITS = re.compile(r"[0-9]")


def _regenerate(tmp_path: Path) -> Path:
    completed = subprocess.run(
        [sys.executable, str(EXAMPLE), "--output", str(tmp_path), "--no-plot"],
        capture_output=True,
        text=True,
        cwd=ROOT,
        timeout=600,
    )
    assert completed.returncode == 0, f"the example failed:\n{completed.stdout}\n{completed.stderr}"
    return tmp_path


@pytest.fixture(scope="module")
def regenerated(tmp_path_factory: pytest.TempPathFactory) -> Path:
    return _regenerate(tmp_path_factory.mktemp("gz_bundle"))


def _compare(left: object, right: object, path: str = "$") -> None:
    if isinstance(right, dict):
        assert isinstance(left, dict), path
        assert set(left) == set(right), f"{path}: keys {sorted(set(left) ^ set(right))}"
        for key in sorted(right):
            _compare(left[key], right[key], f"{path}.{key}")
    elif isinstance(right, list):
        assert isinstance(left, list), path
        assert len(left) == len(right), f"{path}: length {len(left)} != {len(right)}"
        for index, (left_item, right_item) in enumerate(zip(left, right, strict=True)):
            _compare(left_item, right_item, f"{path}[{index}]")
    elif isinstance(right, bool) or right is None:
        assert left is right or left == right, f"{path}: {left!r} != {right!r}"
    elif isinstance(right, (int, float)):
        assert isinstance(left, (int, float)) and not isinstance(left, bool), path
        assert math.isclose(float(left), float(right), rel_tol=1e-9, abs_tol=1e-12), (
            f"{path}: {left!r} != {right!r}"
        )
    else:
        assert left == right, f"{path}: {left!r} != {right!r}"


def _published_results() -> dict:
    return json.loads((BUNDLE / "results.json").read_text())


def test_bundle_is_tracked_with_all_three_artefacts() -> None:
    for name in ("results.json", "report.md", "gz_curves.png"):
        path = BUNDLE / name
        assert path.is_file(), f"{path} is missing; run python examples/gz/run_example.py"
        assert path.stat().st_size > 0


def test_figure_is_a_plot_sized_png() -> None:
    payload = (BUNDLE / "gz_curves.png").read_bytes()
    assert payload[:8] == b"\x89PNG\r\n\x1a\n", "not a PNG"
    width, height = struct.unpack(">II", payload[16:24])
    assert width >= 800 and height >= 300, f"figure is {width}x{height}; the page needs the three panels"
    assert len(payload) > 20_000, "a 3-panel figure of this size should be larger than 20 kB"


def test_version_stamp_matches_the_tree_that_published_it() -> None:
    assert _published_results()["version"] == tensorcad.__version__


def test_published_numbers_reproduce_on_the_current_tree(regenerated: Path) -> None:
    _compare(json.loads((regenerated / "results.json").read_text()), _published_results())


def test_published_page_is_the_page_the_script_writes(regenerated: Path) -> None:
    published = (BUNDLE / "report.md").read_text().splitlines()
    fresh = (regenerated / "report.md").read_text().splitlines()
    assert len(published) == len(fresh), "the report changed shape; regenerate docs/gz/report.md"
    for index, (left, right) in enumerate(zip(published, fresh, strict=True)):
        # Digits are allowed to wobble in the last place between machines (the JSON comparison
        # above is the numeric gate); the words and the table shapes are not.
        assert DIGITS.sub("", left) == DIGITS.sub("", right), f"line {index + 1} changed: {left!r}"


def test_every_claim_in_the_bundle_holds() -> None:
    results = _published_results()
    assert results["checks"]["passed"] is True, results["checks"]["failures"]
    assert not results["checks"]["failures"]
    claims = results["claims"]
    assert len(claims) >= 8, "the page is supposed to carry the whole claim list"
    for entry in claims:
        assert entry["holds"] is True, entry
        assert isinstance(entry["measured"], float)
        assert entry["unit"] in {"m", "relative"}


def test_the_bundle_measures_both_closed_forms_and_a_hull() -> None:
    cases = _published_results()["cases"]
    assert set(cases) == {"box_barge", "log", "hull", "free_surface", "upright"}
    barge, log, hull = cases["box_barge"], cases["log"], cases["hull"]
    # the two closed-form references are sampled at the same heels as the integrated curve
    assert len(barge["gz_m"]) == len(barge["wall_sided_m"]) == len(barge["angles_deg"])
    assert len(log["gz_m"]) == len(log["exact_m"]) == len(log["angles_deg"])
    assert len(hull["gz_m"]) == len(hull["angles_deg"])
    # and the hull curve is a righting-arm curve, not a straight line of numbers
    assert hull["max_gz_m"] > 0.0
    assert hull["angle_of_max_gz_deg"] in hull["angles_deg"]
    if hull["angle_of_max_gz_deg"] == hull["angles_deg"][-1]:
        # a maximum sitting on the last sample is the edge of the range, not a peak, and the
        # range of stability must then say so instead of inventing an angle
        assert hull["range_of_stability_deg"] is None
    assert hull["max_volume_error"] < 1e-10
