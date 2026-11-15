"""Tests for the Verdict oracle (``tensorcad.meshquality.verdict_oracle``).

These tests exercise the *oracle*, not the metric module: they assert that the compiled
Verdict library answers what ``tests/test_meshquality.py`` already asserts internally, so a
regression in either side is visible.  On a host without a C++ build of the oracle the whole
module is skipped -- the run is recorded as "not run", never as a pass.
"""
from __future__ import annotations

import math

import pytest

from tensorcad.meshquality import verdict_oracle
from tensorcad.meshquality.metrics import hex_metrics, tet_metrics, tri_metrics
from tensorcad.meshquality.reference import external_availability

pytestmark = pytest.mark.skipif(
    not verdict_oracle.available(),
    reason="no Verdict oracle built (python tools/verdict_oracle/build.py)",
)

CUBE = [
    (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
    (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 1.0), (0.0, 1.0, 1.0),
]
STRETCHED = [(4.0 * x, y, z) for x, y, z in CUBE]
TRI = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)]
TET = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)]
# A regular tetrahedron: unit edges, so its scaled Jacobian is 1/sqrt(2) in this tree's
# normalisation and 1 in Verdict's.
REGULAR_TET = [
    (0.5, 0.5 * math.sqrt(1.0 / 3.0) * 0.0, 0.0),
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.0),
    (0.5, math.sqrt(3.0) / 2.0, 0.0),
]
REGULAR_TET_3D = [
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.0),
    (0.5, math.sqrt(3.0) / 2.0, 0.0),
    (0.5, math.sqrt(3.0) / 6.0, math.sqrt(2.0 / 3.0)),
]


def test_availability_and_provenance() -> None:
    assert verdict_oracle.available() is True
    assert verdict_oracle.info().startswith("verdict-oracle/")
    provenance = verdict_oracle.provenance()
    assert provenance["built"] is True
    assert provenance["compiler"], "provenance must record the compiler that built the oracle"
    assert provenance["library_sha256"]
    assert external_availability()["verdict_oracle"].startswith("verdict-oracle/")


def test_metric_names_are_stable() -> None:
    for kind in verdict_oracle.KINDS:
        names = verdict_oracle.metric_names(kind)
        assert "scaled_jacobian" in names
        assert "volume" in names or "area" in names
    with pytest.raises(ValueError):
        verdict_oracle.metric_names("wedge")


def test_ideal_element_metrics_match_the_statements_in_the_documentation() -> None:
    quad = verdict_oracle.metrics("quad", CUBE[:4], [(0, 1, 2, 3)])
    assert quad["scaled_jacobian"][0] == pytest.approx(1.0, abs=1e-15)
    assert quad["min_angle"][0] == pytest.approx(90.0, abs=1e-12)
    assert quad["area"][0] == pytest.approx(1.0, abs=1e-15)

    hex_quality = verdict_oracle.metrics("hex", CUBE, [tuple(range(8))])
    assert hex_quality["scaled_jacobian"][0] == pytest.approx(1.0, abs=1e-15)
    assert hex_quality["volume"][0] == pytest.approx(1.0, abs=1e-15)

    # Orthogonal but stretched: shape-only metrics stay 1, the edge ratio sees the stretch.
    stretched = verdict_oracle.metrics("hex", STRETCHED, [tuple(range(8))])
    assert stretched["scaled_jacobian"][0] == pytest.approx(1.0, abs=1e-15)
    assert stretched["edge_ratio"][0] == pytest.approx(4.0, abs=1e-15)

    tri = verdict_oracle.metrics("tri", TRI, [(0, 1, 2)])
    assert tri["area"][0] == pytest.approx(0.5, abs=1e-16)


def test_hex_metrics_agree_with_this_tree_on_the_three_comparable_pairs() -> None:
    cells = [tuple(range(8))]
    for points in (CUBE, STRETCHED):
        ours = hex_metrics(points, cells)
        theirs = verdict_oracle.metrics("hex", points, cells)
        assert ours.metrics["scaled_jacobian"][0] == pytest.approx(
            theirs["scaled_jacobian"][0], abs=1e-14
        )
        assert ours.metrics["edge_ratio"][0] == pytest.approx(theirs["edge_ratio"][0], abs=1e-14)
        assert ours.metrics["volume"][0] == pytest.approx(theirs["volume"][0], abs=1e-14)


def test_triangle_metrics_agree_with_verdict() -> None:
    points = [(0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (0.3, 1.1, 0.0)]
    ours = tri_metrics(points, [(0, 1, 2)])
    theirs = verdict_oracle.metrics("tri", points, [(0, 1, 2)])
    for name in ("scaled_jacobian", "min_angle", "max_angle", "edge_ratio", "area"):
        assert ours.metrics[name][0] == pytest.approx(theirs[name][0], abs=1e-12), name


def test_tetrahedron_scaled_jacobian_differs_by_exactly_one_over_sqrt_two() -> None:
    """The documented normalisation difference, asserted as algebra rather than prose."""
    for points in (TET, REGULAR_TET_3D):
        cells = [(0, 1, 2, 3)]
        ours = tet_metrics(points, cells).metrics["scaled_jacobian"][0]
        theirs = verdict_oracle.metrics("tet", points, cells)["scaled_jacobian"][0]
        assert theirs != 0.0
        assert ours / theirs == pytest.approx(1.0 / math.sqrt(2.0), rel=1e-12)
    # and the volume agrees outright
    ours = tet_metrics(TET, [(0, 1, 2, 3)]).metrics["volume"][0]
    theirs = verdict_oracle.metrics("tet", TET, [(0, 1, 2, 3)])["volume"][0]
    assert ours == pytest.approx(theirs, abs=1e-15)


def test_a_collapsed_hexahedron_is_refused_by_verdict_and_signed_by_this_tree() -> None:
    """Verdict returns its DBL_MAX sentinel; this tree reports a non-positive Jacobian.

    The difference is deliberate (an inverted cell must stay visible), so the test pins both
    behaviours: it fails if this tree starts returning a sentinel, and it fails if Verdict
    starts scoring a collapsed cell as if it were positive.
    """
    collapsed = [(0.0, 0.0, 0.0)] * 8
    ours = hex_metrics(collapsed, [tuple(range(8))]).metrics["scaled_jacobian"][0]
    theirs = verdict_oracle.metrics("hex", collapsed, [tuple(range(8))])["scaled_jacobian"][0]
    assert ours <= 0.0
    assert theirs >= 1e29


def test_unknown_kind_and_bad_node_count_raise() -> None:
    with pytest.raises(ValueError):
        verdict_oracle.metrics("wedge", CUBE, [(0, 1, 2, 3, 4, 5)])
    with pytest.raises(ValueError):
        verdict_oracle.metrics("tri", CUBE, [(0, 1, 2, 3)])