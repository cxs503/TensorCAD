"""Independent-route verification of ``tensorcad.fillet`` (mission m16561).

The module's own tests compare two derivations that both live inside ``tensorcad.fillet``
(Green's-theorem path integral vs the per-corner area closed form).  These tests recompute the
same quantities on routes that share no code with the kernel: a hand-derived rounded-rectangle
closed form, brute-force midpoint integration of the body's indicator function, Gauss-Legendre
quadrature of ``A(r(x))`` along the sweep, and Richardson ratios of the chorded outline.

The example script ``examples/fillet_independent_verification/run_example.py`` is the source of
the published numbers; these tests re-derive them at a smaller grid so a merge gate turns red if
either the kernel or the reference drifts.
"""
from __future__ import annotations

import importlib.util
import math
from pathlib import Path

import pytest

from tensorcad import fillet as F

_SCRIPT = Path(__file__).resolve().parents[1] / "examples" / "fillet_independent_verification" / "run_example.py"
_spec = importlib.util.spec_from_file_location("fillet_independent_verification", _SCRIPT)
assert _spec is not None and _spec.loader is not None
INDEP = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(INDEP)

SQUARE = [(0.0, 0.0), (4.0, 0.0), (4.0, 2.0), (0.0, 2.0)]
L_SECTION = [(0.0, 0.0), (2.0, 0.0), (2.0, 1.0), (1.0, 1.0), (1.0, 2.0), (0.0, 2.0)]


def _centroidal(shape) -> tuple[float, float]:
    moments = shape.moments
    cu, cv = moments.centroid
    return (moments.second_u - moments.area * cu * cu, moments.second_v - moments.area * cv * cv)


@pytest.mark.parametrize("radius", [0.05, 0.2, 0.4, 0.8, 0.999])
def test_rounded_rectangle_matches_a_hand_derived_closed_form(radius: float) -> None:
    """Area, perimeter and both centroidal second moments come from a hand derivation."""
    shape = F.fillet_section(SQUARE, radius)
    ref = INDEP.rounded_rectangle_closed_form(4.0, 2.0, radius)
    assert shape.moments.area == pytest.approx(ref["area"], rel=1e-13)
    assert shape.moments.perimeter == pytest.approx(ref["perimeter"], rel=1e-13)
    i_uu, i_vv = _centroidal(shape)
    assert i_uu == pytest.approx(ref["i_uu"], rel=1e-12)
    assert i_vv == pytest.approx(ref["i_vv"], rel=1e-12)


def test_second_moments_survive_a_brute_force_indicator_integration() -> None:
    """A midpoint rule that only asks "is this point inside?" reproduces I_uu and I_vv."""
    shape = F.fillet_section(SQUARE, 0.4)
    brute = INDEP.brute_force_moments(4.0, 2.0, 0.4, n=400)
    i_uu, i_vv = _centroidal(shape)
    # O(h) on the curved boundary: the brute-force grid, not the kernel, sets this tolerance.
    assert brute["area"] == pytest.approx(shape.moments.area, rel=2e-3)
    assert brute["i_uu"] == pytest.approx(i_uu, rel=2e-3)
    assert brute["i_vv"] == pytest.approx(i_vv, rel=2e-3)


def test_reentrant_section_follows_the_signed_corner_law() -> None:
    base = F.polygon_moments(L_SECTION)
    turns = F.section_turns(L_SECTION)
    for radius in (0.1, 0.3, 0.5):
        shape = F.fillet_section(L_SECTION, radius)
        d_area = sum(
            -math.copysign(1.0, turn) * radius * radius * (math.tan(abs(turn) / 2.0) - abs(turn) / 2.0)
            for turn in turns
            if abs(turn) > 1e-12
        )
        d_perimeter = sum(
            radius * abs(turn) - 2.0 * radius * math.tan(abs(turn) / 2.0)
            for turn in turns
            if abs(turn) > 1e-12
        )
        assert shape.moments.area == pytest.approx(base.area + d_area, rel=1e-13)
        assert shape.moments.perimeter == pytest.approx(base.perimeter + d_perimeter, rel=1e-13)
    # the reflex corner is the only one that adds material
    assert F.fillet_section(L_SECTION, 0.3).moments.area < base.area


def test_the_exact_area_is_the_limit_of_the_chorded_outline() -> None:
    """Richardson over arc_points: the exact value must be the mesh's own limit, order 2."""
    errors = []
    for chords in (4, 8, 16, 32):
        shape = F.fillet_section(SQUARE, 0.4, arc_points=chords)
        sampled = shape.sampled_moments().area
        errors.append((sampled - shape.moments.area) / shape.moments.area)
    ratios = [errors[i] / errors[i + 1] for i in range(len(errors) - 1)]
    assert all(error < 0.0 for error in errors), "chords cut the arc: sampled area is smaller"
    assert sum(ratios) / len(ratios) == pytest.approx(4.0, rel=0.05)


def test_swept_volume_matches_gauss_legendre_quadrature_of_the_area() -> None:
    coefficient = F._area_coefficient(F.section_turns(SQUARE), "fillet")
    base_area = F.polygon_moments(SQUARE).area
    for r0, r1 in ((0.2, 0.6), (0.6, 0.2), (0.1, 0.999), (0.4, 0.4)):
        solid = F.fillet_solid(SQUARE, 3.0, r0, radius_end=r1, arc_points=16, stations=8)
        quad = INDEP.sweep_volume_quadrature(base_area, coefficient, 3.0, r0, r1)
        assert solid.volume_exact == pytest.approx(quad, rel=1e-13)
        assert solid.volume_mesh == pytest.approx(quad, rel=5e-4)
        assert solid.watertight and solid.euler_characteristic == 2
        assert solid.zero_area_faces == 0


def test_circle_limit_is_the_circle_and_the_next_radius_is_refused() -> None:
    square2 = [(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)]
    limit = F.max_fillet_radius(square2)
    assert limit == pytest.approx(1.0, rel=1e-12)
    circle = F.fillet_section(square2, limit)
    assert circle.moments.area == pytest.approx(math.pi, rel=1e-13)
    assert circle.moments.perimeter == pytest.approx(2.0 * math.pi, rel=1e-13)
    with pytest.raises(F.FilletError, match="largest feasible"):
        F.fillet_section(square2, 1.02)


def test_near_limit_sections_stay_simple_and_orientation_independent() -> None:
    cases = {
        "sliver": ([(0.0, 0.0), (4.0, 0.0), (4.0, 0.2), (0.0, 0.2)], 0.0999),
        "rectangle": (SQUARE, 0.999),
        "l_section": (L_SECTION, 0.5),
        "notched": ([(0.0, 0.0), (3.0, 0.0), (3.0, 1.0), (2.0, 1.0),
                     (2.0, 0.4), (1.0, 0.4), (1.0, 1.0), (0.0, 1.0)], 0.29),
    }
    for tag, (points, radius) in cases.items():
        shape = F.fillet_section(points, radius, arc_points=8)
        assert shape.moments.area > 0.0, tag
        assert shape.area_crosscheck_error < 1e-12, tag
        assert shape.sampled_moments().area < shape.moments.area, tag
        clockwise = F.fillet_section(list(reversed(points)), radius, arc_points=8)
        assert clockwise.moments.area == pytest.approx(shape.moments.area, rel=1e-12), tag


def test_translation_conditioning_is_bounded_as_documented() -> None:
    """Absolute-coordinate Green's theorem: error grows like eps * shift / body size."""
    reference = F.fillet_section(SQUARE, 0.4).moments.area
    moved = F.fillet_section([(x + 1e6, y + 1e6) for x, y in SQUARE], 0.4)
    assert abs(moved.moments.area - reference) / reference < 1e-9
    moved_far = F.fillet_section([(x + 1e8, y + 1e8) for x, y in SQUARE], 0.4)
    assert abs(moved_far.moments.area - reference) / reference < 1e-7