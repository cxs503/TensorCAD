"""The cheap gate must not report an ERROR on a tree the expensive gate passes.

The gate has two carriers with different ladders: ``.githooks/pre-commit`` runs
``python -m tensorcad.golden --staged --quick`` (one rung) and CI runs the full ladder.
``inclined_box`` holds its floating position by bisection, so the *level axis is the solve
tolerance*, not a grid: ``coarse`` solves the waterline offset to 3e-5 x scale, ``fine`` to
1e-13 x scale.  Its two invariants declared one hand-picked band each (1e-6 m3 and 1e-8 m),
which the fine solve meets to a ten-thousandth of the band and the coarse solve misses by a
factor of 2600 -- so the hook reported

    ERROR inclined_box invariant_violated identity_volume_equals_requested_displacement:
          residual -3.868e-03 > 1e-06

on every tree the full suite passed.  A gate that is red on green trees is not a gate: the
only way past it was ``TENSORCAD_GOLDEN_SKIP=1`` on every commit that touched geometry,
which is exactly the edit class the hook exists for.

The fix derives both bands from the level's own solve accuracy (:meth:`GoldenCase.invariants`
now receives the acceptance rung) instead of typing one number for two rungs.  These tests
pin the three properties that keep it a gate: the two ladders agree on the ERROR set, the
bands really do track the solve tolerance, and a genuine regression still trips them.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping

import pytest

from tensorcad import golden
from tensorcad.golden_inclined_box import InclinedBoxGolden

ROOT = Path(__file__).resolve().parents[1]
BASELINE = ROOT / "docs" / "golden" / "tensorcad_geometry_golden.json"

#: The band this case declared for both invariants before the level was offered to it.
LEGACY_VOLUME_BAND = 1e-6


def _single_case_baseline(tmp_path: Path, name: str) -> Path:
    """The committed baseline restricted to one case, so a run can register only that case."""
    data = json.loads(BASELINE.read_text(encoding="utf-8"))
    data["cases"] = {name: data["cases"][name]}
    path = tmp_path / "baseline.json"
    path.write_text(json.dumps(data), encoding="utf-8")
    return path


def _error_codes(result: golden.SuiteResult) -> list[tuple[str, str]]:
    return sorted((finding.case, finding.code) for finding in result.errors)


def _invariant_bands(level: str) -> dict[str, tuple[float, float, str]]:
    case = InclinedBoxGolden()
    return case.invariants(case.measure(level), level=level)


class _LevelProbe(golden.GoldenCase):
    """A case that records which rung the engine asserts its invariants on."""

    name = "level_probe"
    title = "records the acceptance level"
    reference_source = "none: the probe measures nothing"

    def __init__(self) -> None:
        self.levels_seen: list[str | None] = []

    def levels(self) -> tuple[str, ...]:
        return ("coarse", "fine")

    def measure(self, level: str) -> dict[str, float]:
        return {"x": 1.0 if level == "coarse" else 1.0 + 1e-6}

    def reference(self) -> dict[str, float]:
        return {"x": 1.0}

    def invariants(
        self, values: Mapping[str, float], *, level: str | None = None
    ) -> dict[str, tuple[float, float, str]]:
        self.levels_seen.append(level)
        return {"identity_x_is_positive": (max(0.0, -values["x"]), 0.0, "x > 0")}


class _LegacyProbe(golden.GoldenCase):
    """A case that still declares invariants with the older one-argument signature."""

    name = "legacy_probe"
    title = "one-argument invariants"
    reference_source = "none: the probe measures nothing"

    def levels(self) -> tuple[str, ...]:
        return ("coarse", "fine")

    def measure(self, level: str) -> dict[str, float]:
        return {"x": 1.0}

    def reference(self) -> dict[str, float]:
        return {"x": 1.0}

    def invariants(self, values: Mapping[str, float]) -> dict[str, tuple[float, float, str]]:  # type: ignore[override]
        return {"identity_x_is_positive": (max(0.0, -values["x"]), 0.0, "x > 0")}


def test_the_quick_hook_and_the_full_ladder_agree_on_the_error_set(tmp_path: Path) -> None:
    """The property the round-1 defect broke: same tree, same verdict, either ladder."""
    path = _single_case_baseline(tmp_path, "inclined_box")
    case = InclinedBoxGolden()
    quick = golden.run_suite(quick=True, baseline_path=path, cases=[case])
    full = golden.run_suite(quick=False, baseline_path=path, cases=[case])

    assert _error_codes(quick) == [], "the hook would block a tree CI passes"
    assert _error_codes(full) == [], "the full ladder must stay clean"
    assert _error_codes(quick) == _error_codes(full)
    assert quick.passed and full.passed


def test_the_hook_really_did_fail_before_the_band_followed_the_level() -> None:
    """The bug was load-bearing, not cosmetic: the coarse residual misses the old band."""
    coarse = _invariant_bands("coarse")
    fine = _invariant_bands("fine")

    residual_coarse = coarse["identity_volume_equals_requested_displacement"][0]
    assert abs(residual_coarse) > LEGACY_VOLUME_BAND, (
        "the legacy 1e-6 m3 band held at the coarse rung, so there was no defect to fix"
    )
    # ...and the derived band holds at that same rung, which is the fix.
    assert abs(residual_coarse) <= coarse["identity_volume_equals_requested_displacement"][1]
    for band in (coarse, fine):
        for key, (residual, tolerance, _basis) in band.items():
            assert abs(residual) <= tolerance, (key, residual, tolerance)


def test_each_band_tracks_the_solve_tolerance_of_its_own_level() -> None:
    """No number is typed for two rungs: the bands move with the level axis."""
    case = InclinedBoxGolden()
    coarse = case.invariants(case.measure("coarse"), level="coarse")
    fine = case.invariants(case.measure("fine"), level="fine")
    ratio = case.SOLVE_TOLERANCE["coarse"] / case.SOLVE_TOLERANCE["fine"]
    assert ratio > 1e6, ratio  # the level axis really is the solve tolerance

    for key in coarse:
        coarse_band = coarse[key][1]
        fine_band = fine[key][1]
        assert fine_band < coarse_band, (key, fine_band, coarse_band)
        # A floor (round-off) is added on top of the solve-accuracy term, so the two bands
        # are not separated by the full solve-tolerance ratio; what must hold is that the
        # band *follows* the tolerance -- loosening one rung's tolerance loosens its band.
        case.SOLVE_TOLERANCE = dict(case.SOLVE_TOLERANCE)
        case.SOLVE_TOLERANCE["coarse"] *= 10.0
        loosened = case.invariants(case.measure("coarse"), level="coarse")
        assert loosened[key][1] == pytest.approx(10.0 * coarse_band, rel=1e-3), (
            key,
            loosened[key][1],
            coarse_band,
        )
        case.SOLVE_TOLERANCE = InclinedBoxGolden.SOLVE_TOLERANCE

    # the basis states where the number came from, so a reader can check it without the code
    basis = coarse["identity_volume_equals_requested_displacement"][2]
    assert "half-width" in basis and "coarse" in basis


def test_the_derived_band_still_trips_on_a_real_regression() -> None:
    """A gate that cannot fail is not a gate: a 2 % displacement error must be caught."""
    case = InclinedBoxGolden()
    values = dict(case.measure("coarse"))
    values["displacement_volume_m3"] *= 1.02
    residual, tolerance, _basis = case.invariants(values, level="coarse")[
        "identity_volume_equals_requested_displacement"
    ]
    assert abs(residual) > tolerance, (residual, tolerance)


def test_a_perturbed_attitude_route_is_caught_through_the_engine(tmp_path: Path) -> None:
    """End to end: a 2 % volume error in the measured module fails a full-ladder run."""
    path = _single_case_baseline(tmp_path, "inclined_box")
    case = InclinedBoxGolden()

    class _Drifted(InclinedBoxGolden):
        def measure(self, level: str) -> dict[str, float]:
            values = super().measure(level)
            values["displacement_volume_m3"] *= 1.02
            return values

    result = golden.run_suite(quick=False, baseline_path=path, cases=[_Drifted()])
    codes = {finding.code for finding in result.errors}
    assert "invariant_violated" in codes, codes
    assert not result.passed
    assert golden.run_suite(quick=False, baseline_path=path, cases=[case]).passed


def test_the_engine_offers_the_acceptance_level_to_the_case(tmp_path: Path) -> None:
    """Quick asserts the first rung, the full ladder its acceptance rung -- not ``None``."""
    path = _single_case_baseline(tmp_path, "inclined_box")
    quick = _LevelProbe()
    golden.run_suite(quick=True, baseline_path=path, baseline_optional=True, cases=[quick])
    assert quick.levels_seen == ["coarse"]

    full = _LevelProbe()
    golden.run_suite(quick=False, baseline_path=path, baseline_optional=True, cases=[full])
    assert full.levels_seen == ["fine"]


def test_a_case_with_the_older_signature_is_still_evaluated(tmp_path: Path) -> None:
    """The adapter must not silently drop a case that never learned about levels."""
    path = _single_case_baseline(tmp_path, "inclined_box")
    result = golden.run_suite(quick=True, baseline_path=path, baseline_optional=True,
                              cases=[_LegacyProbe()])
    probe = result.cases[0]
    invariants = [metric for metric in probe.metrics if metric.kind == "invariant"]
    assert [metric.name for metric in invariants] == ["identity_x_is_positive"]
    assert probe.passed


def test_the_committed_baseline_records_the_level_derived_bands() -> None:
    """The evidence file carries the derived numbers, not the hand-picked ones."""
    data: dict[str, Any] = json.loads(BASELINE.read_text(encoding="utf-8"))
    recorded = data["cases"]["inclined_box"]["invariants"]
    fine = _invariant_bands("fine")
    for key, (residual, tolerance, basis) in fine.items():
        entry = recorded[key]
        assert entry["tolerance"] == pytest.approx(tolerance, rel=1e-12), key
        assert entry["residual"] == pytest.approx(residual, rel=1e-6, abs=1e-15), key
        assert entry["basis"] == basis, key