"""Areas under the exact multi-regime box curve: closed form, kinks, and who they certify.

The tests are aimed at the *integral*, not at the curve: the curve already has its own file.
Inside the wall-sided window the area has a closed-form antiderivative, so the integrator can be
held to it; outside, the only independent statement left is the fundamental theorem (the derived
area must differentiate back to ``GZ``) plus agreement with a different quadrature on a different
curve source -- :mod:`tensorcad.statutory`'s adaptive Simpson on the same callable, and the
heeled-mesh curve it was written for.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from tensorcad import statutory
from tensorcad.box_gz import (
    BoxGzError,
    box_gz,
    box_gz_area,
    box_regime_bounds_deg,
    box_section,
)
from tensorcad.geometry import MeshBuilder
from tensorcad.gz import righting_arm_curve

# (breadth, draft, depth, kg), the same four geometries the curve tests use
GEOMETRIES = {
    "deck first": (2.0, 0.8, 1.2, 0.55),
    "bilge first": (2.0, 0.8, 2.4, 0.6),
    "wide shallow": (4.0, 0.5, 0.8, 0.9),
    "square deep": (3.0, 1.0, 1.5, 1.2),
}

#: Five-point Gauss-Legendre, repeated here on purpose: the test below has to be able to
#: integrate *without* the regime split to show what the split buys.
_NODES = (
    0.0,
    -0.5384693101056831,
    0.5384693101056831,
    -0.9061798459386640,
    0.9061798459386640,
)
_WEIGHTS = (
    0.5688888888888889,
    0.4786286704993665,
    0.4786286704993665,
    0.2369268850561891,
    0.2369268850561891,
)


def wall_sided_area(stop_deg: float, *, breadth: float, draft: float, kg: float) -> float:
    """``GM0 (1 - cos) + BM0/2 (sec + cos - 2)``: the antiderivative of the closed form."""
    bm = breadth**2 / (12.0 * draft)
    gm = 0.5 * draft + bm - kg
    phi = math.radians(stop_deg)
    return gm * (1.0 - math.cos(phi)) + 0.5 * bm * (1.0 / math.cos(phi) + math.cos(phi) - 2.0)


def unsplit_kn_area(geometry: tuple[float, float, float, float], stop_deg: float, panels: int) -> float:
    """Composite Gauss-Legendre over the whole interval -- no split at the kinks."""
    breadth, draft, depth, _ = geometry
    low, high = 0.0, math.radians(stop_deg)
    width = (high - low) / panels
    half = 0.5 * width
    total = 0.0
    for index in range(panels):
        centre = low + (index + 0.5) * width
        for node, weight in zip(_NODES, _WEIGHTS, strict=True):
            total += weight * box_section(
                breadth=breadth, draft=draft, depth=depth, heel_deg=math.degrees(centre + half * node)
            )[2]
    return total * half


def box_mesh(length: float, breadth: float, draft: float, freeboard: float):
    """Closed box, bottom at ``-draft`` and deck at ``+freeboard`` -- the mesh the curve came from."""
    builder = MeshBuilder()
    corner = {}
    for x in (0.0, length):
        for y in (-0.5 * breadth, 0.5 * breadth):
            for z in (-draft, freeboard):
                corner[(x, y, z)] = builder.add_vertex((x, y, z))

    def quad(p0, p1, p2, p3) -> None:
        a, b, c, d = (corner[point] for point in (p0, p1, p2, p3))
        builder.add_triangle(a, b, c)
        builder.add_triangle(a, c, d)

    lo, hi = -0.5 * breadth, 0.5 * breadth
    bottom, top = -draft, freeboard
    left, right = 0.0, length
    quad((left, lo, bottom), (left, hi, bottom), (right, hi, bottom), (right, lo, bottom))
    quad((left, lo, top), (right, lo, top), (right, hi, top), (left, hi, top))
    quad((left, lo, bottom), (right, lo, bottom), (right, lo, top), (left, lo, top))
    quad((left, hi, bottom), (left, hi, top), (right, hi, top), (right, hi, bottom))
    quad((left, lo, bottom), (left, lo, top), (left, hi, top), (left, hi, bottom))
    quad((right, lo, bottom), (right, hi, bottom), (right, hi, top), (right, lo, top))
    return builder.build()


@pytest.mark.parametrize("label", sorted(GEOMETRIES))
@pytest.mark.parametrize("fraction", [0.4, 0.8, 0.95])
def test_inside_the_window_the_area_is_the_closed_form_antiderivative(
    label: str, fraction: float
) -> None:
    """Where the wall-sided form still holds, the new integral must reproduce its antiderivative."""
    breadth, draft, depth, kg = GEOMETRIES[label]
    bounds = box_regime_bounds_deg(breadth=breadth, draft=draft, depth=depth)
    stop_deg = fraction * min(bounds["deck_immersion"], bounds["bilge_emergence"])
    assert stop_deg < min(bounds["deck_immersion"], bounds["bilge_emergence"])
    result = box_gz_area(breadth=breadth, draft=draft, depth=depth, kg=kg, stop_deg=stop_deg)
    expected = wall_sided_area(stop_deg, breadth=breadth, draft=draft, kg=kg)
    assert result.area == pytest.approx(expected, rel=1e-11)
    assert result.quadrature_error < 1e-14 * max(abs(result.kn_area), 1.0)


def test_the_interval_is_cut_where_the_curve_kinks() -> None:
    """The regime bounds are not decoration: they are the split points, and they show up."""
    breadth, draft, depth, kg = GEOMETRIES["deck first"]
    bounds = box_regime_bounds_deg(breadth=breadth, draft=draft, depth=depth)
    whole = box_gz_area(breadth=breadth, draft=draft, depth=depth, kg=kg, stop_deg=89.0)
    # both regime boundaries live below 89 degrees and both walls are lost at the second
    assert whole.breakpoints_deg == pytest.approx(
        (bounds["deck_immersion"], bounds["bilge_emergence"]), rel=0.0
    )
    # a window that straddles nothing inside (30-40) is a single piece
    straddling = box_gz_area(
        breadth=breadth, draft=draft, depth=depth, kg=kg, start_deg=15.0, stop_deg=35.0
    )
    assert straddling.breakpoints_deg == pytest.approx((bounds["deck_immersion"],), rel=0.0)
    inside_one_regime = box_gz_area(
        breadth=breadth, draft=draft, depth=depth, kg=kg, start_deg=25.0, stop_deg=35.0
    )
    assert inside_one_regime.breakpoints_deg == ()
    assert inside_one_regime.panels <= whole.panels


def test_the_split_is_what_makes_the_quadrature_converge() -> None:
    """Without the split the deck-edge kink costs six digits, and more panels do not buy it back."""
    geometry = GEOMETRIES["deck first"]
    breadth, draft, depth, kg = geometry
    reference = box_gz_area(
        breadth=breadth, draft=draft, depth=depth, kg=kg, stop_deg=40.0, tolerance=1e-15
    ).kn_area
    unsplit_256 = unsplit_kn_area(geometry, 40.0, 256)
    # ten times the panels of the split run, and still wrong in the eleventh digit
    assert abs(unsplit_256 - reference) > 1e-12
    split = box_gz_area(breadth=breadth, draft=draft, depth=depth, kg=kg, stop_deg=40.0)
    assert split.heel_evaluations < 5 * 256
    assert abs(split.kn_area - reference) < 1e-15


@pytest.mark.parametrize("stop_deg", [12.0, 30.0, 55.0, 80.0])
def test_the_derived_area_differentiates_back_to_the_curve(stop_deg: float) -> None:
    """The fundamental theorem, as a check that the integral is of *this* curve."""
    breadth, draft, depth, kg = GEOMETRIES["deck first"]
    step = 1e-4
    ahead = box_gz_area(
        breadth=breadth, draft=draft, depth=depth, kg=kg, stop_deg=stop_deg + step
    ).area
    behind = box_gz_area(
        breadth=breadth, draft=draft, depth=depth, kg=kg, stop_deg=stop_deg - step
    ).area
    derivative = (ahead - behind) / (2.0 * math.radians(step))
    gz = box_gz(breadth=breadth, draft=draft, depth=depth, kg=kg, heel_deg=stop_deg).righting_arm
    assert derivative == pytest.approx(gz, abs=1e-9)


@pytest.mark.parametrize("label", sorted(GEOMETRIES))
def test_the_kg_term_is_analytic_and_therefore_exact(label: str) -> None:
    """``KG sin(phi)`` is integrated in closed form, so a KG shift moves the area by an exact amount."""
    breadth, draft, depth, kg = GEOMETRIES[label]
    start, stop, shift = 0.0, 40.0, 0.25
    base = box_gz_area(breadth=breadth, draft=draft, depth=depth, kg=kg, start_deg=start, stop_deg=stop)
    heavier = box_gz_area(
        breadth=breadth, draft=draft, depth=depth, kg=kg + shift, start_deg=start, stop_deg=stop
    )
    expected = shift * (math.cos(math.radians(start)) - math.cos(math.radians(stop)))
    assert heavier.kn_area == pytest.approx(base.kn_area, rel=1e-14)
    assert base.area - heavier.area == pytest.approx(expected, rel=1e-12)
    assert heavier.kg_term - base.kg_term == pytest.approx(expected, rel=1e-12)


@pytest.mark.parametrize("label", sorted(GEOMETRIES))
def test_the_result_is_odd_in_the_interval(label: str) -> None:
    """Heeling the other way is the same integral with the sign flipped -- no hidden branch."""
    breadth, draft, depth, kg = GEOMETRIES[label]
    positive = box_gz_area(breadth=breadth, draft=draft, depth=depth, kg=kg, stop_deg=35.0)
    negative = box_gz_area(
        breadth=breadth, draft=draft, depth=depth, kg=kg, start_deg=-35.0, stop_deg=0.0
    )
    assert negative.area == pytest.approx(-positive.area, rel=1e-12)
    assert negative.breakpoints_deg == pytest.approx(
        tuple(-angle for angle in reversed(positive.breakpoints_deg)), rel=0.0
    )


def test_the_area_is_additive_across_the_breakpoints() -> None:
    """Splitting a window at an arbitrary heel must not change the sum."""
    breadth, draft, depth, kg = GEOMETRIES["deck first"]
    whole = box_gz_area(breadth=breadth, draft=draft, depth=depth, kg=kg, stop_deg=45.0).area
    parts = sum(
        box_gz_area(breadth=breadth, draft=draft, depth=depth, kg=kg, start_deg=low, stop_deg=high).area
        for low, high in ((0.0, 21.801409486351808), (21.801409486351808, 42.0), (42.0, 45.0))
    )
    assert parts == pytest.approx(whole, rel=1e-13)


@pytest.mark.parametrize(
    ("start_deg", "stop_deg"),
    [(0.0, 30.0), (0.0, 40.0), (30.0, 40.0), (0.0, 60.0), (0.0, 89.0), (20.0, 25.0)],
)
def test_statutory_agrees_at_large_heel_where_no_antiderivative_exists(
    start_deg: float, stop_deg: float
) -> None:
    """Two independent quadratures on two independent curve sources, out to 89 degrees.

    The IS Code areas are what :mod:`tensorcad.statutory` computes for a rule check; the only
    reference it had stopped at the wall-sided limit.  This is that reference extended.
    """
    breadth, draft, depth, kg = GEOMETRIES["deck first"]
    exact = box_gz_area(
        breadth=breadth, draft=draft, depth=depth, kg=kg, start_deg=start_deg, stop_deg=stop_deg
    )
    callable_curve = lambda angle: box_gz(  # noqa: E731 - the point is that it is a plain callable
        breadth=breadth, draft=draft, depth=depth, kg=kg, heel_deg=angle
    ).righting_arm
    theirs = statutory.area_under_gz_curve(callable_curve, start_deg, stop_deg)
    assert theirs == pytest.approx(exact.area, rel=1e-12)
    assert exact.quadrature_error <= 1e-12 * max(abs(exact.kn_area), 1.0)


def test_the_heeled_mesh_curve_agrees_too() -> None:
    """The same areas off the mesh integrator the reference was built to certify."""
    breadth, draft, depth, kg = GEOMETRIES["deck first"]
    mesh = box_mesh(10.0, breadth, draft, depth - draft)
    angles = [float(angle) for angle in range(0, 86)]
    curve = righting_arm_curve(mesh, kg, angles_deg=angles)
    table = (list(curve.angles_deg), list(curve.gz))
    for start_deg, stop_deg in ((0.0, 30.0), (0.0, 40.0), (0.0, 85.0)):
        exact = box_gz_area(
            breadth=breadth, draft=draft, depth=depth, kg=kg, start_deg=start_deg, stop_deg=stop_deg
        ).area
        theirs = statutory.area_under_gz_curve(table, start_deg, stop_deg)
        assert theirs == pytest.approx(exact, rel=1e-4)


def test_the_reported_error_bounds_the_observed_one() -> None:
    """A tighter tolerance buys nothing, and the estimate never understates the damage by much."""
    breadth, draft, depth, kg = GEOMETRIES["deck first"]
    loose = box_gz_area(breadth=breadth, draft=draft, depth=depth, kg=kg, stop_deg=40.0)
    tight = box_gz_area(
        breadth=breadth, draft=draft, depth=depth, kg=kg, stop_deg=40.0, tolerance=1e-16
    )
    assert tight.panels >= loose.panels
    observed = abs(tight.area - loose.area)
    assert observed <= max(loose.quadrature_error * 1.0e3, 1e-14 * max(abs(loose.area), 1.0))


def test_a_default_call_is_cheap_enough_to_loop_over_a_hull() -> None:
    """The IS Code wants three areas and the weather criterion wants four; this must be usable."""
    breadth, draft, depth, kg = GEOMETRIES["deck first"]
    result = box_gz_area(breadth=breadth, draft=draft, depth=depth, kg=kg, stop_deg=40.0)
    assert result.heel_evaluations <= 200
    assert result.quadrature_error < 1e-15


@pytest.mark.parametrize(
    "kwargs",
    [
        {"stop_deg": 0.0},
        {"start_deg": 30.0, "stop_deg": 30.0},
        {"start_deg": 40.0, "stop_deg": 30.0},
        {"stop_deg": 89.5},
        {"stop_deg": 90.0},
        {"start_deg": -90.0},
        {"stop_deg": math.nan},
        {"stop_deg": math.inf},
        {"kg": math.nan},
        {"tolerance": 0.0},
        {"tolerance": -1e-9},
    ],
)
def test_the_model_refuses_what_it_cannot_represent(kwargs: dict) -> None:
    breadth, draft, depth, kg = GEOMETRIES["deck first"]
    arguments = {"breadth": breadth, "draft": draft, "depth": depth, "kg": kg, **kwargs}
    with pytest.raises(BoxGzError):
        box_gz_area(**arguments)


def test_a_draft_deeper_than_the_depth_is_refused() -> None:
    with pytest.raises(BoxGzError):
        box_gz_area(breadth=2.0, draft=1.5, depth=1.2, kg=0.6)

ARTIFACT = Path(__file__).resolve().parents[1] / "artifacts" / "m4873" / "box_gz_areas.json"
PAGE = Path(__file__).resolve().parents[1] / "docs" / "box_gz_areas.md"


def test_the_committed_evidence_still_makes_its_claims() -> None:
    """Re-apply the harness thresholds to the committed json, and pin the numbers to the page."""
    payload = json.loads(ARTIFACT.read_text())
    thresholds = payload["thresholds"]
    page = PAGE.read_text()
    results = payload["results"]
    assert results
    assert (
        max(row["relative_difference"] for result in results for row in result["antiderivative"])
        <= thresholds["antiderivative"]
    )
    assert (
        max(row["relative_difference"] for result in results for row in result["areas"])
        <= thresholds["statutory"]
    )
    assert (
        max(row["relative_difference"] for result in results for row in result["mesh_areas"])
        <= thresholds["mesh"]
    )
    assert min(
        result["controls"]["wall_sided_past_limit_relative_difference"] for result in results
    ) >= thresholds["control"]
    # every regulated area in the artifact is on the page, so a re-run without a doc sync shows
    for result in results:
        for row in result["areas"]:
            # the residuals table renders the exact area at twelve decimals
            assert f"{row['exact']:.12f}" in page


def test_the_artifact_is_what_the_code_computes_today() -> None:
    """A stale artifact is a lie with a timestamp; recompute the headline number."""
    payload = json.loads(ARTIFACT.read_text())
    result = payload["results"][0]
    spec = result["spec"]
    headline = box_gz_area(
        breadth=spec["breadth"],
        draft=spec["draft"],
        depth=result["depth"],
        kg=spec["kg"],
        start_deg=0.0,
        stop_deg=40.0,
    ).area
    committed = next(row for row in result["areas"] if row["window_deg"] == [0.0, 40.0])["exact"]
    assert headline == pytest.approx(committed, rel=1e-13)
