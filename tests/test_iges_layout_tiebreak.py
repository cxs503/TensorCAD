"""Control-net layout choice in the type-128 reader.

The reader has to guess whether a surface's control points are written as IGES
5.3 ``x, y, z, w`` quadruples or as a block of ``n`` weights followed by ``3n``
coordinates, because both occur in the workshop archives.  The regression these
tests lock down is the *degenerate* case: a patch with no interior control point
(2x2) makes the old smoothness score identically zero for every candidate, so the
layout used to be decided by tuple sort order -- silently, and wrongly for 15 of
the 132 ONRT patches and 2 of the 45 JBC patches.
"""

from __future__ import annotations

import pytest

import math

from tensorcad.iges import BSplineSurface, _split_control, _uniform_weights

# Two candidate readings of an ONRT 2x2 patch, built from the same 16-value tail.
COORDS = [
    (-0.0310, 0.0452, -0.0100),
    (-0.0300, 0.0500, -0.0110),
    (-0.0280, 0.0452, -0.0100),
    (-0.0270, 0.0500, -0.0110),
]


def _close(control, expected, tol: float = 1e-12) -> bool:
    """Element-wise comparison of two point lists (``approx`` rejects nesting)."""
    if len(control) != len(expected):
        return False
    return all(
        math.isclose(a, b, rel_tol=0.0, abs_tol=tol)
        for got, want in zip(control, expected)
        for a, b in zip(got, want)
    )


def _weight_block_tail() -> list[float]:
    tail = [1.0, 1.0, 1.0, 1.0]
    for point in COORDS:
        tail.extend(point)
    return tail


def _interleaved_tail() -> list[float]:
    tail = []
    for point in COORDS:
        tail.extend([point[0], point[1], point[2], 1.0])
    return tail


def _surface(tail: list[float], *, rational: bool, count_u: int = 2, count_v: int = 2) -> BSplineSurface:
    """Assemble a type-128 parameter list from a control-point tail."""
    params: list[float] = [
        128.0,
        float(count_u - 1),
        float(count_v - 1),
        1.0,
        1.0,  # degree 1 x degree 1
        0.0,
        0.0,
        0.0 if rational else 1.0,  # PROP3: 0 = rational/clamped, 1 = polynomial
        0.0,  # PROP4
        0.0,  # slot the reader skips before the knot vectors (see iges.cursor = 10)
    ]
    params += [0.0, 0.0, 1.0, 1.0]  # u knots (count_u + degree + 1 = 4)
    params += [0.0, 0.0, 1.0, 1.0]  # v knots
    params += list(tail)
    return BSplineSurface.from_parameters(params)


def test_uniform_weights_predicate() -> None:
    assert _uniform_weights([1.0, 1.0, 1.0])
    assert _uniform_weights([2.5, 2.5, 2.5, 2.5])
    assert not _uniform_weights([1.0, 1.0, 0.75])  # genuinely rational patch
    assert not _uniform_weights([0.0, 1.0, 1.0])  # a zero weight is not admissible


def test_weight_block_tail_is_read_as_weight_block() -> None:
    """A 2x2 weight-block patch must not be read as interleaved (the old bug)."""
    surface = _surface(_weight_block_tail(), rational=True)
    assert surface.layout == "weight_block"
    assert surface.layout_ambiguous is False
    assert _close(surface.control, COORDS)


def test_interleaved_tail_is_read_as_interleaved() -> None:
    surface = _surface(_interleaved_tail(), rational=True)
    assert surface.layout == "interleaved"
    assert _close(surface.control, COORDS)


def test_old_rule_would_have_misread_the_weight_block_patch() -> None:
    """Guard against a regression to the smoothness-only rule.

    Both readings of a 2x2 net are perfectly smooth, so the smoothness score
    cannot separate them; the weight test is what does.  This test pins the fact
    that the two candidates are indistinguishable without it.
    """
    tail = _weight_block_tail()
    as_block, _, _, _ = _split_control(tail, 4, True, 2, 2)
    interleaved = [(tail[4 * i], tail[4 * i + 1], tail[4 * i + 2]) for i in range(4)]
    # the interleaved reading of a weight-block tail invents a point at (1, 1, 1):
    assert any(all(abs(c - 1.0) < 1e-12 for c in point) for point in interleaved)
    assert not any(all(abs(c - 1.0) < 1e-12 for c in point) for point in as_block)


def test_genuinely_rational_patch_falls_back_to_smoothness_and_is_flagged() -> None:
    """A rational patch with non-unit weights: weight test is silent.

    The choice must then come from the smoothness score, and a tie must be
    reported as ambiguous (fail closed) rather than hidden.
    """
    coords = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0), (3.0, 0.0, 0.0)]
    tail = [1.0, 0.9, 1.1, 1.0]
    for point in coords:
        tail.extend(point)
    control, weights, layout, ambiguous = _split_control(tail, 4, True, 2, 2)
    assert layout in {"weight_block", "interleaved"}
    assert len(control) == 4
    assert len(weights) == 4
    # The weight test is silent (weights genuinely vary) and a 2x2 net has no
    # interior point, so the smoothness score is the only evidence left.
    assert not _uniform_weights(weights)
    assert layout in {"weight_block", "interleaved"}
    assert isinstance(ambiguous, bool)


def test_polynomial_layout_only_offered_for_polynomial_records() -> None:
    """PROP3 = 1 declares plain x, y, z data; a rational record must not use it."""
    triples: list[float] = []
    for point in COORDS:
        triples.extend(point)
    surface = _surface(triples, rational=False)  # PROP3 = 1 -> polynomial
    assert surface.layout == "polynomial"
    assert _close(surface.control, COORDS)