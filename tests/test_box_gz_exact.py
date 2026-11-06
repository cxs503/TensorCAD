"""The exact multi-regime box solution: what it must reproduce, and what must refuse.

The tests are aimed at the *reference*, not at the integrators: a reference that is wrong
everywhere would still agree with itself.  So the checks here are (a) the wall-sided closed form
inside its window, (b) the classical bilge-emergence heel ``atan(2T/B)`` on a geometry where the
deck cannot interfere, (c) continuity where the section changes regime, (d) odd symmetry in heel,
and (e) volume conservation.  Only then is it used, in
``tools/verify_box_large_heel.py``, to certify the two mesh integrators.
"""

from __future__ import annotations

import math

import pytest

from tensorcad.box_gz import (
    BoxGzError,
    box_gz,
    box_gz_curve,
    box_regime_bounds_deg,
    box_section,
    box_vanishing_angle_deg,
)
from tensorcad.geometry import MeshBuilder
from tensorcad.gz import righting_arm_curve, wall_sided_righting_arm

# (breadth, draft, depth, kg) -- labelled by which regime boundary comes first
GEOMETRIES = {
    "deck first": (2.0, 0.8, 1.2, 0.55),
    "bilge first": (2.0, 0.8, 2.4, 0.6),
    "wide shallow": (4.0, 0.5, 0.8, 0.9),
    "square deep": (3.0, 1.0, 1.5, 1.2),
}


def box_mesh(length: float, breadth: float, draft: float, freeboard: float):
    """Closed box with its bottom at ``-draft`` and its deck at ``+freeboard``."""
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


# ---------------------------------------------------------------------------
# the reference against the two closed forms it generalises
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("label", sorted(GEOMETRIES))
def test_inside_its_window_the_exact_solution_is_the_wall_sided_form(label: str) -> None:
    """A generalisation that disagrees with the special case it extends is simply wrong."""
    breadth, draft, depth, kg = GEOMETRIES[label]
    limit = min(
        box_regime_bounds_deg(breadth=breadth, draft=draft, depth=depth)["deck_immersion"],
        box_regime_bounds_deg(breadth=breadth, draft=draft, depth=depth)["bilge_emergence"],
    )
    bm = breadth**2 / (12.0 * draft)
    gm = 0.5 * draft + bm - kg
    for heel in (0.0, 1.0, 0.25 * limit, 0.5 * limit, 0.75 * limit, 0.99 * limit):
        exact = box_gz(
            breadth=breadth, draft=draft, depth=depth, kg=kg, heel_deg=heel
        ).righting_arm
        closed = wall_sided_righting_arm(heel, gm=gm, bm=bm)
        assert exact == pytest.approx(closed, abs=1e-13), (label, heel)


@pytest.mark.parametrize("label", sorted(GEOMETRIES))
def test_the_small_angle_slope_is_the_metacentric_height(label: str) -> None:
    """``GZ'(0) = GM`` -- one scalar that ties the curve to `tensorcad.stability`."""
    breadth, draft, depth, kg = GEOMETRIES[label]
    bm = breadth**2 / (12.0 * draft)
    gm = 0.5 * draft + bm - kg
    heel = 1e-3
    exact = box_gz(breadth=breadth, draft=draft, depth=depth, kg=kg, heel_deg=heel).righting_arm
    assert exact / math.sin(math.radians(heel)) == pytest.approx(gm, rel=1e-8)


def test_bilge_emergence_matches_the_classical_formula_when_the_deck_keeps_out() -> None:
    """``tan(phi) = 2T/B`` is the textbook bilge-emergence heel; the solver must find it.

    The tall barge has enough freeboard that the deck edge is nowhere near the waterline at
    38.66 degrees, so this is the classical configuration and the two numbers must agree.
    """
    breadth, draft, depth, _ = GEOMETRIES["bilge first"]
    bounds = box_regime_bounds_deg(breadth=breadth, draft=draft, depth=depth)
    classical = math.degrees(math.atan(2.0 * draft / breadth))
    assert bounds["bilge_emergence"] == pytest.approx(classical, abs=1e-9)
    assert bounds["deck_immersion"] > bounds["bilge_emergence"]
    assert box_gz(
        breadth=breadth, draft=draft, depth=depth, kg=0.6, heel_deg=classical + 1.0
    ).regime == "bilge_emerged"


def test_a_deck_that_immerses_first_delays_the_bilge_emergence() -> None:
    """Not a tautology: capping volume above the deck raises the waterline and keeps the bilge wet."""
    breadth, draft, depth, _ = GEOMETRIES["deck first"]
    bounds = box_regime_bounds_deg(breadth=breadth, draft=draft, depth=depth)
    classical = math.degrees(math.atan(2.0 * draft / breadth))
    assert bounds["deck_immersion"] < bounds["bilge_emergence"]
    assert bounds["bilge_emergence"] > classical + 1.0


# ---------------------------------------------------------------------------
# the reference's own consistency
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("label", sorted(GEOMETRIES))
@pytest.mark.parametrize("heel", [0.0, 5.0, 21.8, 22.0, 38.66, 40.0, 60.0, 85.0])
def test_the_section_always_holds_the_upright_displacement(label: str, heel: float) -> None:
    breadth, draft, depth, _ = GEOMETRIES[label]
    _, area, _, _ = box_section(breadth=breadth, draft=draft, depth=depth, heel_deg=heel)
    assert area == pytest.approx(breadth * draft, rel=1e-14)


@pytest.mark.parametrize("label", sorted(GEOMETRIES))
def test_the_curve_is_continuous_across_every_regime_boundary(label: str) -> None:
    """A mis-assembled regime shows up as a jump here, not as a 1 % error somewhere."""
    breadth, draft, depth, kg = GEOMETRIES[label]
    bounds = box_regime_bounds_deg(breadth=breadth, draft=draft, depth=depth)
    for name in ("deck_immersion", "bilge_emergence", "both_walls_lost"):
        boundary = bounds[name]
        if not math.isfinite(boundary) or boundary <= 0.0:
            continue
        below = box_gz(
            breadth=breadth, draft=draft, depth=depth, kg=kg, heel_deg=boundary - 1e-6
        ).righting_arm
        above = box_gz(
            breadth=breadth, draft=draft, depth=depth, kg=kg, heel_deg=boundary + 1e-6
        ).righting_arm
        # a mis-assembled regime jumps by centimetres here; the true curve moves ~1e-8 m
        # over a 2e-6 degree stride
        assert abs(above - below) < 1e-7, (label, name, below, above)


@pytest.mark.parametrize("label", sorted(GEOMETRIES))
@pytest.mark.parametrize("heel", [3.0, 25.0, 47.0, 80.0])
def test_the_curve_is_odd_in_heel(label: str, heel: float) -> None:
    """Starboard and port are mirror images, so ``GZ(-phi) = -GZ(phi)`` exactly."""
    breadth, draft, depth, kg = GEOMETRIES[label]
    starboard = box_gz(breadth=breadth, draft=draft, depth=depth, kg=kg, heel_deg=heel)
    port = box_gz(breadth=breadth, draft=draft, depth=depth, kg=kg, heel_deg=-heel)
    assert port.righting_arm == pytest.approx(-starboard.righting_arm, abs=1e-12)
    assert port.regime == starboard.regime


@pytest.mark.parametrize("label", sorted(GEOMETRIES))
def test_the_regimes_are_entered_in_the_order_the_bounds_promise(label: str) -> None:
    breadth, draft, depth, kg = GEOMETRIES[label]
    bounds = box_regime_bounds_deg(breadth=breadth, draft=draft, depth=depth)
    first, second = sorted((bounds["deck_immersion"], bounds["bilge_emergence"]))
    early = box_gz(
        breadth=breadth, draft=draft, depth=depth, kg=kg, heel_deg=0.5 * first
    ).regime
    middle = box_gz(
        breadth=breadth, draft=draft, depth=depth, kg=kg, heel_deg=0.5 * (first + second)
    ).regime
    late = box_gz(breadth=breadth, draft=draft, depth=depth, kg=kg, heel_deg=85.0).regime
    assert early == "wall_sided"
    assert middle in {"deck_immersed", "bilge_emerged"}
    assert late == "bilge_and_deck"


def test_the_vanishing_angle_is_where_the_curve_crosses_zero() -> None:
    breadth, draft, depth, kg = GEOMETRIES["wide shallow"]
    vanish = box_vanishing_angle_deg(breadth=breadth, draft=draft, depth=depth, kg=kg)
    assert 0.0 < vanish < 89.0
    assert box_gz(
        breadth=breadth, draft=draft, depth=depth, kg=kg, heel_deg=vanish
    ).righting_arm == pytest.approx(0.0, abs=1e-9)
    assert box_gz(
        breadth=breadth, draft=draft, depth=depth, kg=kg, heel_deg=vanish - 1.0
    ).righting_arm > 0.0
    assert box_gz(
        breadth=breadth, draft=draft, depth=depth, kg=kg, heel_deg=vanish + 1.0
    ).righting_arm < 0.0


def test_a_stable_barge_can_have_no_vanishing_angle_below_89_degrees() -> None:
    breadth, draft, depth, kg = GEOMETRIES["deck first"]
    assert math.isinf(
        box_vanishing_angle_deg(breadth=breadth, draft=draft, depth=depth, kg=kg)
    )


def test_an_unstable_barge_vanishes_immediately() -> None:
    """GM < 0 means the curve starts downwards; that is angle zero, reported as such."""
    breadth, draft, depth, _ = GEOMETRIES["deck first"]
    assert box_vanishing_angle_deg(breadth=breadth, draft=draft, depth=depth, kg=2.0) == 0.0


def test_the_peak_of_the_curve_is_interior_for_a_soft_barge() -> None:
    breadth, draft, depth, kg = GEOMETRIES["wide shallow"]
    points = box_gz_curve(breadth=breadth, draft=draft, depth=depth, kg=kg)
    arms = [point.righting_arm for point in points]
    peak = arms.index(max(arms))
    assert 0 < peak < len(arms) - 1


# ---------------------------------------------------------------------------
# refusals
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "kwargs",
    [
        {"breadth": 0.0, "draft": 0.8, "depth": 1.2, "kg": 0.5, "heel_deg": 10.0},
        {"breadth": 2.0, "draft": -0.8, "depth": 1.2, "kg": 0.5, "heel_deg": 10.0},
        {"breadth": 2.0, "draft": 0.8, "depth": 0.8, "kg": 0.5, "heel_deg": 10.0},
        {"breadth": 2.0, "draft": 0.8, "depth": 1.2, "kg": float("nan"), "heel_deg": 10.0},
        {"breadth": 2.0, "draft": 0.8, "depth": 1.2, "kg": 0.5, "heel_deg": 90.0},
        {"breadth": 2.0, "draft": 0.8, "depth": 1.2, "kg": 0.5, "heel_deg": float("inf")},
    ],
)
def test_invalid_inputs_are_refused(kwargs: dict) -> None:
    with pytest.raises(BoxGzError):
        box_gz(**kwargs)


def test_an_empty_angle_list_is_refused() -> None:
    with pytest.raises(BoxGzError):
        box_gz_curve(breadth=2.0, draft=0.8, depth=1.2, kg=0.5, angles_deg=())


def test_a_draft_deeper_than_the_depth_is_refused() -> None:
    with pytest.raises(BoxGzError):
        box_regime_bounds_deg(breadth=2.0, draft=1.4, depth=1.2)


# ---------------------------------------------------------------------------
# then, and only then, the integrators
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("heel", [5.0, 25.0, 40.0, 60.0, 80.0])
def test_the_newton_integrator_matches_the_exact_solution_at_large_heel(heel: float) -> None:
    """``tensorcad.gz`` had no independent reference past the wall-sided window until now."""
    breadth, draft, depth, kg = GEOMETRIES["deck first"]
    mesh = box_mesh(10.0, breadth, draft, depth - draft)
    exact = box_gz(
        breadth=breadth, draft=draft, depth=depth, kg=kg, heel_deg=heel
    ).righting_arm
    measured = righting_arm_curve(mesh, kg, angles_deg=(heel,)).gz[0]
    assert measured == pytest.approx(exact, abs=1e-10)


@pytest.mark.parametrize("heel", [40.0, 60.0, 80.0])
def test_the_bisection_integrator_matches_the_exact_solution_at_large_heel(heel: float) -> None:
    """Same for the bracket-and-bisect route carried in from PR #40."""
    from tensorcad import bisection_righting_arm_curve

    breadth, draft, depth, kg = GEOMETRIES["deck first"]
    mesh = box_mesh(10.0, breadth, draft, depth - draft)
    exact = box_gz(
        breadth=breadth, draft=draft, depth=depth, kg=kg, heel_deg=heel
    ).righting_arm
    curve = bisection_righting_arm_curve(mesh, kg=kg, heels_deg=(heel - 0.5, heel))
    measured = curve.points[-1].righting_arm
    assert measured == pytest.approx(exact, abs=1e-9)


def test_a_millimetre_of_kg_is_visible() -> None:
    """Negative control: the comparison above must be able to fail."""
    breadth, draft, depth, kg = GEOMETRIES["deck first"]
    base = box_gz(breadth=breadth, draft=draft, depth=depth, kg=kg, heel_deg=40.0).righting_arm
    shifted = box_gz(
        breadth=breadth, draft=draft, depth=depth, kg=kg + 1e-3, heel_deg=40.0
    ).righting_arm
    assert abs(base - shifted) == pytest.approx(1e-3 * math.sin(math.radians(40.0)), rel=1e-6)
    assert abs(base - shifted) > 1e-4

def _evidence() -> dict:
    import json
    from pathlib import Path

    path = Path(__file__).resolve().parents[1] / "artifacts" / "m4864" / "box_large_heel.json"
    if not path.exists():
        pytest.skip(f"the recorded evidence {path.name} is not in this checkout")
    return json.loads(path.read_text())


def test_the_page_tables_are_the_ones_the_evidence_yields() -> None:
    """The doc's result tables must equal what the recorded evidence regenerates.

    ``tools/verify_box_large_heel.py`` writes those four table blocks into
    ``docs/box_large_heel.md`` in place; this reads the committed JSON back through the same
    formatter, so an edit to either the page or the evidence shows up here rather than as a
    figure on a page nobody re-ran.
    """
    import importlib.util
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    page = (root / "docs" / "box_large_heel.md").read_text()
    spec = importlib.util.spec_from_file_location(
        "verify_box_large_heel", root / "tools" / "verify_box_large_heel.py"
    )
    try:
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    except ImportError as error:  # pragma: no cover - matplotlib is an optional extra
        pytest.skip(f"the verification tool needs matplotlib: {error}")
    blocks = [
        block
        for block in module.markdown_tables(_evidence()["specs"]).strip("\n").split("\n\n")
        if block.strip()
    ]
    assert len(blocks) == 4
    for block in blocks:
        assert block in page, block.splitlines()[0]


def test_the_evidence_keeps_the_claims_the_changelog_makes() -> None:
    """The quoted agreement figures, re-checked against the JSON they were read off."""
    payload = _evidence()
    results = payload["specs"]
    assert payload["mission"] == "m4864"
    assert len(results) == 4
    seen = {regime for result in results for regime in result["regimes_present"]}
    assert seen == {"wall_sided", "deck_immersed", "bilge_emerged", "bilge_and_deck"}
    assert payload["worst"]["newton"] < 1e-9
    assert payload["worst"]["bisect"] < 1e-9
    assert payload["worst"]["wall_sided_window"] < 1e-9
    for result in results:
        assert result["volume_residual_relative"] < 1e-12
        # negative controls: each failure mode is visible well above those tolerances
        assert result["control_kg_1mm_shift_m"] > 1e-4
        assert result["control_displacement_half_percent_effect_m"] > 1e-4
    # the one geometry whose deck never interferes must reproduce atan(2T/B) exactly
    deck_free = [r for r in results if "bilge_emerged" in r["regimes_present"]]
    assert len(deck_free) == 1
    assert abs(deck_free[0]["bilge_minus_classical_deg"]) < 1e-4
