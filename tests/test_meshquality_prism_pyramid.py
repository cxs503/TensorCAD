"""Prisms (wedge6) and pyramids (5 nodes): closed forms, second routes, and the gate.

The two shapes are the ones a ship CFD mesh reaches for and the library did not measure:
prism layers on a hull surface and the pyramids a hex/tet transition leaves behind. What
is pinned here, in the same order as ``tests/test_meshquality.py``:

* closed-form values on the ideal prism and the ideal pyramid (the *definitions* are
  right, including what "ideal" can mean for a pyramid at all);
* agreement with the independent route in :mod:`tensorcad.meshquality.reference`
  (finite-difference Jacobians, Cayley-Menger volumes, ``atan2`` angles) on meshes that
  are *not* affine, so a Jacobian that is constant by accident cannot pass;
* the behaviour the gate needs: a clean boundary layer passes, a collapsed prism fails
  and names the element, and a kind the mesh does not contain is vacuous rather than an
  error.
"""
from __future__ import annotations

import math

import pytest

from tensorcad.meshquality import (
    DEFAULT_GATE,
    NODES_PER_CELL,
    measure,
    prism_metrics,
    pyramid_metrics,
)
from tensorcad.meshquality.io import load_spec, spec_from_mapping, write_spec
from tensorcad.meshquality.reference import reference_prisms, reference_pyramids

IDEAL_PRISM = (
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.0),
    (0.5, math.sqrt(3.0) / 2.0, 0.0),
    (0.0, 0.0, 1.0),
    (1.0, 0.0, 1.0),
    (0.5, math.sqrt(3.0) / 2.0, 1.0),
)
IDEAL_PYRAMID = (
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.0),
    (1.0, 1.0, 0.0),
    (0.0, 1.0, 0.0),
    (0.5, 0.5, 1.0),
)
#: A right prism whose top face is slid in x: not affine, so the Jacobian varies inside it.
SHEARED_PRISM = (
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.0),
    (0.0, 1.0, 0.0),
    (0.6, 0.0, 1.0),
    (1.6, 0.0, 1.0),
    (0.6, 1.0, 1.0),
)
#: A right prism squeezed towards a wedge: the boundary-layer collapse case.
WEDGE_PRISM = (
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.0),
    (0.0, 1.0, 0.0),
    (0.5, 0.0, 1.0),
    (1.0, 0.0, 1.0),
    (0.0, 1.0, 1.0),
)
#: A pyramid with a non-planar base (the fourth node lifted out of the base plane).
WARPED_PYRAMID = (
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.0),
    (1.0, 1.0, 0.0),
    (0.0, 1.0, 0.3),
    (0.5, 0.5, 1.0),
)


def _single(metrics, points, cell):
    return {name: values[0] for name, values in metrics(list(points), [cell]).metrics.items()}


def _prism(points):
    return _single(prism_metrics, points, (0, 1, 2, 3, 4, 5))


def _pyramid(points):
    return _single(pyramid_metrics, points, (0, 1, 2, 3, 4))


# --------------------------------------------------------------------------------------
# closed forms: the definitions
# --------------------------------------------------------------------------------------
def test_ideal_prism_scores_the_closed_form_values() -> None:
    q = _prism(IDEAL_PRISM)
    assert q["volume"] == pytest.approx(math.sqrt(3.0) / 4.0, rel=0, abs=1e-15)
    assert q["min_sicn"] == pytest.approx(1.0, abs=1e-12)
    # scaled Jacobian of the ideal prism is sin(60 deg): the edges of one corner are not
    # orthogonal, and unlike the SICN unit this metric is not renormalised.
    assert q["scaled_jacobian"] == pytest.approx(math.sqrt(3.0) / 2.0, abs=1e-12)
    assert q["min_angle"] == pytest.approx(60.0, abs=1e-10)
    assert q["max_angle"] == pytest.approx(90.0, abs=1e-10)
    assert q["skewness"] == pytest.approx(0.0, abs=1e-12)
    assert q["edge_ratio"] == pytest.approx(1.0, abs=1e-12)


def test_ideal_pyramid_scores_the_closed_form_values() -> None:
    q = _pyramid(IDEAL_PYRAMID)
    assert q["volume"] == pytest.approx(1.0 / 3.0, abs=1e-15)
    assert q["min_sicn"] == pytest.approx(1.0, abs=1e-12)
    assert q["scaled_jacobian"] == pytest.approx(math.sqrt(2.0 / 3.0), abs=1e-12)
    # The apex angle of a side face: 2 asin(0.5 / sqrt(1.5)).
    assert q["min_angle"] == pytest.approx(
        math.degrees(2.0 * math.asin(0.5 / math.sqrt(1.5))), abs=1e-10
    )
    assert q["edge_ratio"] == pytest.approx(math.sqrt(1.5), abs=1e-12)


def test_equiangle_skewness_has_no_ideal_square_pyramid_at_unit_height() -> None:
    """A square pyramid cannot have 60 deg side-face angles at height 1; at 1/sqrt2 it can.

    This is why ``DEFAULT_LIMITS['pyramid']`` does not gate skewness: the "ideal element"
    of the equiangle definition is a shape, not the shape a mesh generator builds.
    """
    assert _pyramid(IDEAL_PYRAMID)["skewness"] == pytest.approx(0.19683858159631004, abs=1e-12)
    h = 1.0 / math.sqrt(2.0)
    flat = (
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (1.0, 1.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.5, 0.5, h),
    )
    q = _pyramid(flat)
    assert q["skewness"] == pytest.approx(0.0, abs=1e-12)
    assert q["min_angle"] == pytest.approx(60.0, abs=1e-9)
    assert q["max_angle"] == pytest.approx(90.0, abs=1e-9)
    assert q["volume"] == pytest.approx(h / 3.0, abs=1e-15)


def test_aspect_is_visible_in_sicn_and_not_in_scaled_jacobian() -> None:
    """A stretched prism keeps its corner angles but loses conditioning (hex behaviour)."""
    stretched = (
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, 0.0, 3.0),
        (1.0, 0.0, 3.0),
        (0.0, 1.0, 3.0),
    )
    q = _prism(stretched)
    assert q["scaled_jacobian"] == pytest.approx(1.0, abs=1e-12)
    assert q["min_sicn"] < 0.7
    assert q["volume"] == pytest.approx(1.5, abs=1e-15)


def test_pyramid_apex_is_excluded_because_it_is_singular() -> None:
    """Every pyramid map has a singular apex; the metric set must not pretend otherwise."""
    q = _pyramid(IDEAL_PYRAMID)
    assert q["min_sicn"] == pytest.approx(1.0, abs=1e-12)
    # Integrating the Jacobian at the apex would give exactly zero for every pyramid.
    from tensorcad.meshquality.metrics import _det3, _pyramid_jacobian

    apex = _det3(_pyramid_jacobian([tuple(map(float, p)) for p in IDEAL_PYRAMID], 0.5, 0.5, 1.0))
    assert apex == pytest.approx(0.0, abs=1e-15)


# --------------------------------------------------------------------------------------
# invariance
# --------------------------------------------------------------------------------------
def _rotate(points, axis_angle, shift, scale=1.0):
    """Rotate about z by ``axis_angle``, scale, then translate (rigid + similarity)."""
    c, s = math.cos(axis_angle), math.sin(axis_angle)
    out = []
    for x, y, z in points:
        out.append(
            (
                scale * (c * x - s * y) + shift[0],
                scale * (s * x + c * y) + shift[1],
                scale * z + shift[2],
            )
        )
    return tuple(out)


@pytest.mark.parametrize("kind", ["prism", "pyramid"])
def test_metrics_are_invariant_under_rigid_motion_and_scaling(kind: str) -> None:
    points = IDEAL_PRISM if kind == "prism" else IDEAL_PYRAMID
    moved = _rotate(points, 0.7, (3.0, -2.0, 5.0), scale=4.5)
    base = _prism(points) if kind == "prism" else _pyramid(points)
    other = _prism(moved) if kind == "prism" else _pyramid(moved)
    for name, value in base.items():
        # lengths scale, volumes scale with the cube, the ratios and Jacobian metrics do not
        factor = {"volume": 4.5**3, "min_edge": 4.5, "max_edge": 4.5}.get(name, 1.0)
        assert other[name] == pytest.approx(value * factor, rel=1e-11), name


def test_inverted_pyramid_is_reported_negative_not_clamped() -> None:
    inverted = tuple((x, y, -z) for x, y, z in IDEAL_PYRAMID)
    q = _pyramid(inverted)
    assert q["min_sicn"] < 0.0
    assert q["volume"] > 0.0  # the magnitude is the volume, the sign is the orientation


def test_flat_prism_does_not_raise() -> None:
    flat = (
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
    )
    q = _prism(flat)
    assert q["volume"] == 0.0
    assert q["min_sicn"] == pytest.approx(0.0, abs=1e-12)


# --------------------------------------------------------------------------------------
# the independent route
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize(
    "points,kind",
    [
        (IDEAL_PRISM, "prism"),
        (SHEARED_PRISM, "prism"),
        (WEDGE_PRISM, "prism"),
        (IDEAL_PYRAMID, "pyramid"),
        (WARPED_PYRAMID, "pyramid"),
    ],
)
def test_library_agrees_with_the_independent_route(points, kind: str) -> None:
    library = _prism(points) if kind == "prism" else _pyramid(points)
    if kind == "prism":
        reference = reference_prisms(list(points), [(0, 1, 2, 3, 4, 5)])
    else:
        reference = reference_pyramids(list(points), [(0, 1, 2, 3, 4)])
    for name, values in reference.items():
        tolerance = 1e-9 if name in {"scaled_jacobian", "min_sicn"} else 1e-12
        assert library[name] == pytest.approx(values[0], abs=tolerance), name


def test_sheared_prism_is_not_affine_so_the_route_is_not_vacuous() -> None:
    """The Jacobian really varies inside the sheared prism; a constant one would pass."""
    from tensorcad.meshquality.metrics import _det3, _prism_jacobian

    pts = [tuple(map(float, p)) for p in SHEARED_PRISM]
    bottom = _det3(_prism_jacobian(pts, 0.0, 0.0, 0.0))
    top = _det3(_prism_jacobian(pts, 0.0, 0.0, 1.0))
    assert bottom == pytest.approx(1.0, abs=1e-12)
    assert top == pytest.approx(1.0, abs=1e-12)
    middle = _det3(_prism_jacobian(pts, 0.0, 0.0, 0.5))
    assert middle == pytest.approx(1.0, abs=1e-12)
    # The *shape* varies even though the volume element does not: the shear makes the
    # column lengths differ, which is what scaled_jacobian and min_sicn see.
    assert _prism(SHEARED_PRISM)["min_sicn"] < 1.0
    assert _prism(IDEAL_PRISM)["min_sicn"] == pytest.approx(1.0, abs=1e-12)


# --------------------------------------------------------------------------------------
# mesh level: measure(), the exchange format, and the gate
# --------------------------------------------------------------------------------------
def _prism_layer_x(count: int = 3):
    """A stack of ``count`` prisms above a unit-square footprint (one column wide)."""
    points: list[tuple[float, float, float]] = []
    prisms: list[tuple[int, ...]] = []
    for layer in range(count):
        z0, z1 = float(layer), float(layer + 1)
        base = len(points)
        points.extend(
            [
                (0.0, 0.0, z0),
                (1.0, 0.0, z0),
                (0.0, 1.0, z0),
                (0.0, 0.0, z1),
                (1.0, 0.0, z1),
                (0.0, 1.0, z1),
            ]
        )
        prisms.append((base, base + 1, base + 2, base + 3, base + 4, base + 5))
    return tuple(points), prisms


def test_measure_handles_a_prism_layer_and_reports_one_element_per_cell() -> None:
    points, prisms = _prism_layer_x(3)
    quality = measure(points, {"prism": prisms})
    assert len(quality["prism"]) == 3
    assert quality["prism"]["volume"] == pytest.approx((0.5, 0.5, 0.5), abs=1e-15)
    payload = quality.to_dict()
    assert payload["prism"]["count"] == 3
    assert payload["prism"]["summary"]["min_sicn"]["count"] == 3


def test_node_table_is_the_single_source_of_truth() -> None:
    assert NODES_PER_CELL["prism"] == 6
    assert NODES_PER_CELL["pyramid"] == 5
    with pytest.raises(ValueError, match="prism cells need 6 nodes"):
        measure(IDEAL_PRISM, {"prism": [(0, 1, 2, 3, 4)]})
    with pytest.raises(ValueError, match="pyramid cells need 5 nodes"):
        measure(IDEAL_PYRAMID, {"pyramid": [(0, 1, 2, 3, 4, 5)]})
    with pytest.raises(ValueError, match="pyramid cell .* indexes outside"):
        measure(IDEAL_PYRAMID, {"pyramid": [(0, 1, 2, 3, 9)]})


def test_exchange_format_round_trips_prisms_and_pyramids(tmp_path) -> None:
    points = IDEAL_PRISM + IDEAL_PYRAMID
    payload = {
        "points": [list(p) for p in points],
        "elements": {"prism": [[0, 1, 2, 3, 4, 5]], "pyramid": [[6, 7, 8, 9, 10]]},
        "metadata": {"source": "test"},
    }
    spec = spec_from_mapping(payload)
    path = write_spec(tmp_path / "mixed.json", spec)
    reloaded = load_spec(path)
    assert reloaded.elements["prism"] == ((0, 1, 2, 3, 4, 5),)
    assert reloaded.elements["pyramid"] == ((6, 7, 8, 9, 10),)
    with pytest.raises(ValueError, match="prism cell 0 has 5 nodes"):
        spec_from_mapping(
            {"points": [list(p) for p in points], "elements": {"prism": [[0, 1, 2, 3, 4]]}}
        )


def test_gate_passes_a_clean_boundary_layer_and_fails_a_collapsed_prism() -> None:
    points, prisms = _prism_layer_x(3)
    clean = measure(points, {"prism": prisms})
    assert DEFAULT_GATE.evaluate(clean).passed

    collapsed_points = list(points)
    # Collapse the top triangle of the last prism onto an edge (zero area, bad angles):
    # nodes 15,16,17 are that prism's top face.
    collapsed_points[15] = (0.0, 0.0, 3.0)
    collapsed_points[16] = (1.0, 0.0, 3.0)
    collapsed_points[17] = (0.5, 0.0, 3.0)
    collapsed = measure(collapsed_points, {"prism": prisms})
    report = DEFAULT_GATE.evaluate(collapsed)
    assert not report.passed
    failing = {check.metric for check in report.failures}
    assert {"min_sicn", "min_angle"} <= failing
    assert all(check.kind == "prism" for check in report.failures)


def test_gate_is_vacuous_for_a_kind_the_mesh_does_not_contain() -> None:
    points, prisms = _prism_layer_x(1)
    clean = measure(points, {"prism": prisms})
    report = DEFAULT_GATE.evaluate(clean)
    assert report.passed
    assert all(check.kind == "prism" for check in report.checks)


def test_cli_gates_a_prism_mesh_and_returns_the_spec_verdict(tmp_path, capsys) -> None:
    from tensorcad.meshquality.cli import main

    points, prisms = _prism_layer_x(2)
    spec = {
        "points": [list(p) for p in points],
        "elements": {"prism": [list(cell) for cell in prisms]},
    }
    path = tmp_path / "layer.json"
    write_spec(path, spec_from_mapping(spec))
    assert main(["--input", str(path), "--json", str(tmp_path / "out.json"), "--quiet"]) == 0
    assert (tmp_path / "out.json").exists()

    broken = [list(p) for p in points]
    broken[3] = broken[4] = broken[5] = [0.0, 0.0, 1.0]  # zero-thickness first prism
    broken_spec = {"points": broken, "elements": {"prism": [list(cell) for cell in prisms]}}
    broken_path = tmp_path / "broken.json"
    write_spec(broken_path, spec_from_mapping(broken_spec))
    assert main(["--input", str(broken_path), "--quiet"]) == 1
    assert "check(s) failed" in capsys.readouterr().err