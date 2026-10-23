"""Tests for the parameter validation layer."""
from __future__ import annotations

import pytest

from tensorcad.parameters import (
    FLOAT,
    INT,
    ParameterSpec,
    RatioSpec,
    check_ratios,
    require,
    validate,
)

SPECS = (
    ParameterSpec("length", FLOAT, minimum=0.0, maximum=100.0, exclusive_minimum=True),
    ParameterSpec("nx", INT, minimum=3.0, maximum=100.0),
)
RATIOS = (RatioSpec("l_over_two", "length", "half", minimum=1.0, maximum=100.0, guidance="stub"),)


@pytest.mark.parametrize("value", [1.0, 1, 1e-9, 100.0])
def test_accepts_valid_floats(value: float) -> None:
    assert SPECS[0].check(value) is None


@pytest.mark.parametrize("value", [0.0, -1.0, float("nan"), float("inf"), "3", None, True, [1]])
def test_rejects_bad_floats(value: object) -> None:
    assert SPECS[0].check(value) is not None


@pytest.mark.parametrize("value", [3, 4, 100])
def test_accepts_valid_ints(value: int) -> None:
    assert SPECS[1].check(value) is None


@pytest.mark.parametrize("value", [3.0, 4.5, True, "4", 2, 101, float("nan")])
def test_rejects_bad_ints_including_integral_floats(value: object) -> None:
    """A float count is a division that went wrong upstream, not a rounding issue."""
    assert SPECS[1].check(value) is not None


def test_messages_quote_the_unit_and_the_offending_value() -> None:
    message = ParameterSpec("nx", INT, minimum=3.0, unit="stations").check(2)
    assert message == "nx must be >= 3 stations, got 2"
    assert "must be > 0 m" in ParameterSpec("length", FLOAT, minimum=0.0, exclusive_minimum=True).check(0.0)


def test_validate_collects_every_problem_and_flags_unknowns() -> None:
    errors = validate(SPECS, {"length": -1.0, "nx": 1})
    assert len(errors) == 2
    assert "missing parameter 'nx'" in validate(SPECS, {"length": 1.0})[0]
    assert "unknown parameters: beam, zz" in validate(SPECS, {"length": 1.0, "nx": 4, "beam": 1, "zz": 2})[0]


def test_require_raises_one_combined_error() -> None:
    with pytest.raises(ValueError, match=r"length.*; .*nx"):
        require(SPECS, {"length": 0.0, "nx": 2})
    require(SPECS, {"length": 5.0, "nx": 5})


def test_unknown_kind_is_a_programming_error_not_a_data_error() -> None:
    with pytest.raises(ValueError, match="unknown kind"):
        ParameterSpec("x", "complex").check(1.0)


def test_ratio_checks_are_advisory() -> None:
    assert check_ratios(RATIOS, {"length": 10.0, "half": 5.0}) == []
    problem = check_ratios(RATIOS, {"length": 10.0, "half": 0.01})
    assert problem and "l_over_two" in problem[0] and "outside" in problem[0]
    # A bad denominator is the bounds check's business, not the ratio's.
    assert check_ratios(RATIOS, {"length": 10.0, "half": 0.0}) == []
    assert check_ratios(RATIOS, {"length": 10.0, "half": float("inf")}) == []


def test_bounds_stop_an_accidental_oom_before_allocation() -> None:
    from tensorcad import HullParameters, generate_hull

    with pytest.raises(ValueError, match="nx must be <="):
        generate_hull(HullParameters(nx=4002))
    with pytest.raises(ValueError, match="nz must be <="):
        generate_hull(HullParameters(nz=514))
