"""SUBOFF appendages (AFF-8): the sail, the cruciform fins, and the full model.

What is pinned here is external, not self-referential:

* the section values computed by TensorLBM's transcription of the DARPA report
  (``suboff_cad._sail_half_thickness_np`` / ``_naca4_thickness_np``, read
  2026-10-24) — so a change to the polynomials is a change against the source;
* containment labels returned by the *actual* torch predicates
  (``suboff_sail_contains_points`` / ``suboff_fins_contain_points``) at a set of
  hand-picked points, including the tight ones near the cap and the fin's
  leading edge;
* a verbatim transcription of those predicates used as an oracle over 20,000
  random points, so the comparison is not limited to the hand-picked set;
* volumes from quadrature over the published cross-sections, independent of the
  implicit field and of the sampler.

The one thing the tests cannot pin is the *magnitude* of the implicit function
(the sail's cap and the fin's thickness surface have no closed-form nearest
point), so they measure it instead: the gap along the analytic surface normal.
"""

from __future__ import annotations

import math

import pytest

from tensorcad.geometry_field import (
    IMPLICIT_FUNCTION,
    FieldError,
    GeometryField,
    SampledField,
    to_mesh,
)
from tensorcad.suboff import SUBOFF_LENGTH_M, suboff_radius_fraction_array, suboff_reference_geometry
from tensorcad.suboff_appendages import (
    _DECK_FT,
    _FIN_R_INNER_FT,
    _FIN_TRAILING_FT,
    _LENGTH_FT,
    _SAIL_BODY_TOP_FT,
    _SAIL_HALF_WIDTH_FT,
    _SAIL_X_CENTER_FT,
    _SAIL_X_END_FT,
    _SAIL_X_START_FT,
    SuboffAppendagesField,
    SuboffFinsField,
    SuboffFullField,
    SuboffSailField,
    _fin_chord_coefficients,
    _fin_planform_distance_ft,
    _naca_half_thickness_array,
    _sail_half_width_array,
    _sail_half_width_slope_array,
    suboff_fin_chord,
    suboff_sail_half_width,
)

#: Sail half width (ft) at seven stations, from TensorLBM ``_sail_half_thickness_np``.
TENSORLBM_SAIL_HALF_WIDTH = {
    3.05: 0.04904483522377447,
    3.2: 0.10644549991291676,
    3.4: 0.109375,
    3.55: 0.109375,
    3.8: 0.0994715017078891,
    4.0: 0.06692398822021756,
    4.2: 0.014061932204276293,
}

#: NACA half thickness (ft) at seven chord positions, from TensorLBM ``_naca4_thickness_np``.
TENSORLBM_NACA_HALF_THICKNESS = {
    0.02: 0.03932962554685719,
    0.1: 0.07804677373039919,
    0.3: 0.1000287773232838,
    0.5: 0.08815875333428597,
    0.7: 0.060653511877967026,
    0.9: 0.022816421191197506,
    0.98: 0.004799850348000401,
}

#: Containment labels from the torch predicates at the reference length, for points
#: given in feet in the DARPA frame. ``(x, y, z, sail_inside, fins_inside)``.
TENSORLBM_POINTS_FT = [
    (3.6, 0.0, 1.0, True, False),  # sail body, on the centreplane
    (3.6, 0.05, 0.5, True, False),  # sail body, off centre
    (3.6, 0.2, 0.5, False, False),  # beside the sail (|y| > a)
    (3.6, 0.0, 1.53, True, False),  # under the elliptical cap
    (3.6, 0.0, 1.60, False, False),  # above the cap (top is 1.5625 ft)
    (3.6, 0.05, 1.53, True, False),  # cap, off centre
    (3.6, 0.05, 1.60, False, False),  # above the cap, off centre
    (3.05, 0.0, 1.0, True, False),  # entrance segment, centreline
    (4.2, 0.0, 1.0, True, False),  # exit segment, centreline
    (3.0, 0.0, 1.0, False, False),  # forward of the sail
    (12.9, 0.3, 0.0, False, True),  # side fin, mid-blade
    (12.9, 0.3, 0.05, False, True),  # side fin, inside the thickness
    (12.9, 0.3, 0.3, False, False),  # side fin, outside the thickness
    (12.9, 0.05, 0.3, False, True),  # vertical fin, near its root
    (13.0, 0.0, 0.4, False, True),  # vertical fin, on the centreplane
    (12.0, 0.4, 0.0, False, False),  # forward of that span's leading edge
    (12.5, 0.6, 0.0, False, False),  # 3.1e-2 ft forward of the leading edge
    (13.2, 0.6, 0.0, False, False),  # aft of the trailing edge
]

#: Volume of the full model (m^3) from quadrature over the published sections, with
#: the parts buried inside the hull removed. Computed by the example script; pinned
#: here so the sampler is judged against the formulas rather than against itself.
FULL_MODEL_VOLUME_M3 = 0.70754062


def _feet(length: float = SUBOFF_LENGTH_M) -> float:
    return _LENGTH_FT / length


def _to_metres(points_ft: list[tuple[float, float, float]], length: float = SUBOFF_LENGTH_M) -> list:
    feet = _feet(length)
    return [(x / feet, y / feet, z / feet) for x, y, z in points_ft]


# --------------------------------------------------------------------------- #
# The published sections
# --------------------------------------------------------------------------- #


def test_the_sail_section_matches_its_source() -> None:
    feet = _feet()
    for x_ft, expected in TENSORLBM_SAIL_HALF_WIDTH.items():
        assert float(_sail_half_width_array(x_ft)) == pytest.approx(expected, abs=1e-12)
        assert suboff_sail_half_width(x_ft / feet) == pytest.approx(expected / feet, rel=1e-15)


def test_the_sail_section_closes_at_both_tips() -> None:
    assert float(_sail_half_width_array(_SAIL_X_START_FT)) == 0.0
    assert float(_sail_half_width_array(_SAIL_X_END_FT)) == 0.0
    assert float(_sail_half_width_array(0.5 * _SAIL_X_START_FT)) == 0.0
    assert float(_sail_half_width_array(_SAIL_X_END_FT + 0.1)) == 0.0


def test_the_sail_segments_join_without_a_step() -> None:
    """C0 at every joint, at the 1e-9 ft level the polynomials can hold."""
    for joint in (3.358507, 3.559028):
        left = float(_sail_half_width_array(joint - 1e-9))
        right = float(_sail_half_width_array(joint + 1e-9))
        assert abs(left - right) < 1e-7
        assert right == pytest.approx(_SAIL_HALF_WIDTH_FT, rel=1e-6)


def test_the_slope_is_the_derivative_of_the_section() -> None:
    grid = [3.04, 3.15, 3.30, 3.5, 3.7, 3.9, 4.15, 4.23]
    for x_ft in grid:
        step = 1e-7
        numeric = (
            float(_sail_half_width_array(x_ft + step)) - float(_sail_half_width_array(x_ft - step))
        ) / (2.0 * step)
        assert float(_sail_half_width_slope_array(x_ft)) == pytest.approx(numeric, rel=1e-5, abs=1e-5)


def test_the_fin_section_matches_its_source() -> None:
    feet = _feet()
    for s, expected in TENSORLBM_NACA_HALF_THICKNESS.items():
        assert float(_naca_half_thickness_array(s)) == pytest.approx(expected, abs=1e-15)
    # The chord law at the root and tip radii, and the sweep that produces it.
    root = suboff_fin_chord(_FIN_R_INNER_FT / feet)
    assert root == pytest.approx((-0.466308 * _FIN_R_INNER_FT + 0.88859) / feet, rel=1e-15)
    assert root == pytest.approx(0.26018, abs=1e-5)
    tip = suboff_fin_chord(0.825 / feet)
    assert tip == pytest.approx(0.15358, abs=1e-5)


def test_the_planform_distance_is_the_distance_to_the_quadrilateral() -> None:
    """The ``(x, r)`` planform is a trapezoid; its distance must be exact, not a gap."""
    intercept, r_outer = _fin_chord_coefficients(1.0)
    corners = [
        (_FIN_TRAILING_FT - (-0.466308 * _FIN_R_INNER_FT + intercept), _FIN_R_INNER_FT),
        (_FIN_TRAILING_FT, _FIN_R_INNER_FT),
        (_FIN_TRAILING_FT, r_outer),
        (_FIN_TRAILING_FT - (-0.466308 * r_outer + intercept), r_outer),
    ]
    boundary = []
    for index, (ax, ar) in enumerate(corners):
        bx, br = corners[(index + 1) % len(corners)]
        for step in range(501):
            t = step / 500.0
            boundary.append((ax + t * (bx - ax), ar + t * (br - ar)))
    worst = 0.0
    for x_ft in (12.2, 12.4, 12.6, 12.9, 13.1, 12.0):
        for r_ft in (0.0, 0.05, 0.1, 0.3, 0.5, 0.7, 0.85):
            brute = min(math.dist((x_ft, r_ft), point) for point in boundary)
            value = _fin_planform_distance_ft(x_ft, r_ft, 1.0)
            assert abs(value) == pytest.approx(brute, abs=1e-3)
            worst = max(worst, abs(abs(value) - brute))
    assert worst < 1e-3


# --------------------------------------------------------------------------- #
# Containment: the published predicate
# --------------------------------------------------------------------------- #


def _tensorlbm_sail_predicate(points, length: float, scale: float) -> list[bool]:
    """Verbatim transcription of TensorLBM's ``suboff_sail_contains_points``."""
    feet = _LENGTH_FT / length
    inverse = 1.0 / scale
    out = []
    for x_m, y_m, z_m in points:
        x = _SAIL_X_CENTER_FT + (x_m * feet - _SAIL_X_CENTER_FT) * inverse
        y = y_m * feet * inverse
        z = _DECK_FT + (z_m * feet - _DECK_FT) * inverse
        zmax, y_tmp = _SAIL_HALF_WIDTH_FT, _SAIL_BODY_TOP_FT
        inside = False
        if 3.032986 < x < 3.358507:  # m1
            d1 = 3.0720 * (x - 3.032986)
            c1 = 1.0 - (d1 - 1.0) ** 4 * (4.0 * d1 + 1.0)
            b1 = (1.0 / 3.0) * d1 * d1 * (d1 - 1.0) ** 3
            a1 = 2.0 * d1 * (d1 - 1.0) ** 4
            half1 = zmax * math.sqrt(max(2.094759 * a1 + 0.2071781 * b1 + c1, 0.0))
            body = 0.0 < z <= y_tmp and -half1 < y < half1
            cap = y_tmp < z < y_tmp + half1 / 2.0 and abs(y) < math.sqrt(
                max(half1 * half1 - (2.0 * (z - y_tmp)) ** 2, 0.0)
            )
            inside = body or cap
        if 3.358507 < x <= 3.559028:  # m2
            body = 0.0 < z <= y_tmp and -zmax < y < zmax
            cap = y_tmp < z < y_tmp + zmax / 2.0 and abs(y) < math.sqrt(
                max(zmax**2 - (2.0 * (z - y_tmp)) ** 2, 0.0)
            )
            inside = body or cap
        if 3.559028 < x <= 4.241319:  # m3
            e3 = (4.241319 - x) / 0.6822917
            f3 = e3 - 1.0
            half3 = zmax * (
                2.238361 * e3 * f3**4 + 3.106529 * e3 * e3 * f3**3 + 1.0 - f3**4 * (4.0 * e3 + 1.0)
            )
            body = 0.0 < z <= y_tmp and -half3 < y < half3
            cap = y_tmp < z < y_tmp + half3 / 2.0 and abs(y) < math.sqrt(
                max(half3 * half3 - (2.0 * (z - y_tmp)) ** 2, 0.0)
            )
            inside = body or cap
        out.append(inside)
    return out


def _tensorlbm_fins_predicate(points, length: float, scale: float) -> list[bool]:
    """Verbatim transcription of TensorLBM's ``suboff_fins_contain_points``."""
    feet = _LENGTH_FT / length
    intercept, r_outer = _fin_chord_coefficients(scale)
    out = []
    for x_m, y_m, z_m in points:
        x, y, z = x_m * feet, abs(y_m * feet), abs(z_m * feet)
        chord_y = -0.466308 * y + intercept
        chord_z = -0.466308 * z + intercept
        s_y = (x - _FIN_TRAILING_FT) / chord_y + 1.0
        s_z = (x - _FIN_TRAILING_FT) / chord_z + 1.0

        def half(s: float) -> float:
            return scale * float(_naca_half_thickness_array(min(max(s, 0.0), 1.0)))

        fins_y = _FIN_R_INNER_FT < y < r_outer and z < half(s_y) and 0.0 < s_y < 1.0
        fins_z = _FIN_R_INNER_FT < z < r_outer and y < half(s_z) and 0.0 < s_z < 1.0
        out.append(fins_y or fins_z)
    return out


def test_containment_matches_the_labels_the_torch_predicates_returned() -> None:
    feet = _feet()
    points = _to_metres([(x, y, z) for x, y, z, _, _ in TENSORLBM_POINTS_FT])
    sail = SuboffSailField()
    fins = SuboffFinsField()
    assert sail.contains(points) == [row[3] for row in TENSORLBM_POINTS_FT]
    assert fins.contains(points) == [row[4] for row in TENSORLBM_POINTS_FT]
    assert len(points) == 18
    # ... and the field's sign agrees with containment, so a sampler sees the same solid.
    for value, label in zip(sail.signed_distance(points), sail.contains(points), strict=True):
        assert (value < 0.0) == label
    for value, label in zip(fins.signed_distance(points), fins.contains(points), strict=True):
        assert (value < 0.0) == label
    # The metric round trip really is the DARPA frame.
    assert points[0][0] == pytest.approx(3.6 / feet)


def test_containment_matches_the_transcribed_predicate_everywhere() -> None:
    """20,000 random points, three scales, against the verbatim oracle."""
    numpy = pytest.importorskip("numpy")
    generator = numpy.random.default_rng(20261024)
    count = 20_000
    points = numpy.column_stack(
        [
            generator.uniform(0.70, 4.45, count),
            generator.uniform(-0.45, 0.45, count),
            generator.uniform(-0.05, 0.60, count),
        ]
    ).tolist()
    for scale in (1.0, 0.85, 1.2):
        sail = SuboffSailField(scale=scale)
        fins = SuboffFinsField(scale=scale)
        assert sail.contains(points) == _tensorlbm_sail_predicate(points, sail.length, scale)
        assert fins.contains(points) == _tensorlbm_fins_predicate(points, fins.length, scale)
        appendages = SuboffAppendagesField(sail_scale=scale, fin_scale=scale)
        oracle = [
            a or b
            for a, b in zip(
                _tensorlbm_sail_predicate(points, sail.length, scale),
                _tensorlbm_fins_predicate(points, fins.length, scale),
                strict=True,
            )
        ]
        assert appendages.contains(points) == oracle


# --------------------------------------------------------------------------- #
# The magnitude, measured rather than claimed
# --------------------------------------------------------------------------- #


def _sail_surface_samples() -> tuple[list, list]:
    """Analytic sail surface points (metres) and their outward unit normals."""
    feet = _feet()
    points: list[tuple[float, float, float]] = []
    normals: list[tuple[float, float, float]] = []
    for step in range(1, 40):
        x_ft = _SAIL_X_START_FT + (_SAIL_X_END_FT - _SAIL_X_START_FT) * step / 40.0
        width = float(_sail_half_width_array(x_ft))
        slope = float(_sail_half_width_slope_array(x_ft))
        for side in (1.0, -1.0):
            for fraction in (0.15, 0.5, 0.85):
                z_ft = _SAIL_BODY_TOP_FT * fraction
                norm = math.hypot(slope, 1.0)
                points.append((x_ft / feet, side * width / feet, z_ft / feet))
                normals.append((-slope / norm, side / norm, 0.0))
            for index in range(1, 12):
                angle = math.pi * index / 12.0
                y_ft = side * width * math.cos(angle)
                z_ft = _SAIL_BODY_TOP_FT + 0.5 * width * math.sin(angle)
                # rho = |y| below the body top and sqrt(y^2 + (2 rise)^2) above it,
                # so on the cap the gradient is (-a', cos, 2 sin) in unit width.
                norm = math.hypot(slope, math.cos(angle), 2.0 * math.sin(angle))
                points.append((x_ft / feet, y_ft / feet, z_ft / feet))
                normals.append((-slope / norm, side * math.cos(angle) / norm, 2.0 * math.sin(angle) / norm))
    return points, normals


def _fin_surface_samples() -> tuple[list, list]:
    """Analytic fin surface points (metres) and their outward unit normals."""
    feet = _feet()
    intercept, r_outer = _fin_chord_coefficients(1.0)
    points: list[tuple[float, float, float]] = []
    normals: list[tuple[float, float, float]] = []
    for span_step in range(10):
        radius = _FIN_R_INNER_FT + (r_outer - _FIN_R_INNER_FT) * (0.1 + 0.8 * span_step / 9.0)
        chord = -0.466308 * radius + intercept
        for chord_step in range(14):
            s = 0.1 + 0.8 * chord_step / 13.0
            x_ft = _FIN_TRAILING_FT + (s - 1.0) * chord
            half = float(_naca_half_thickness_array(s))
            slope = float(
                _naca_half_thickness_array(min(s + 1e-6, 1.0))
                - _naca_half_thickness_array(max(s - 1e-6, 0.0))
            ) / (min(s + 1e-6, 1.0) - max(s - 1e-6, 0.0))
            ds_dx = 1.0 / chord
            ds_dr = -(x_ft - _FIN_TRAILING_FT) * -0.466308 / (chord * chord)
            for sign in (1.0, -1.0):
                # blade spanning in y (thickness in z), then in z (thickness in y)
                for axis in ("y", "z"):
                    if axis == "y":
                        point = (x_ft, radius, sign * half)
                        normal = (-slope * ds_dx, -slope * ds_dr, sign)
                    else:
                        point = (x_ft, sign * half, radius)
                        normal = (-slope * ds_dx, sign, -slope * ds_dr)
                    norm = math.sqrt(sum(component * component for component in normal))
                    points.append((point[0] / feet, point[1] / feet, point[2] / feet))
                    normals.append(tuple(component / norm for component in normal))
    return points, normals


@pytest.mark.parametrize("offset", [1e-4, 5e-4, 2e-3])
def test_the_gap_along_the_surface_normal_is_first_order_accurate(offset: float) -> None:
    """``f(p + t n) = t`` is what ``f = g / |grad g|`` is *for*; measure how close."""
    sail_points, sail_normals = _sail_surface_samples()
    fin_points, fin_normals = _fin_surface_samples()
    sail = SuboffSailField()
    fins = SuboffFinsField()

    def worst(points: list, normals: list, field: GeometryField, tolerance: float) -> float:
        probes = [
            (p[0] + offset * n[0], p[1] + offset * n[1], p[2] + offset * n[2])
            for p, n in zip(points, normals, strict=True)
        ]
        values = field.signed_distance(probes)
        errors = [abs(value - offset) for value in values]
        median = sorted(errors)[len(errors) // 2]
        assert median < tolerance, f"median {median:.3e} m at offset {offset:.1e} m"
        return max(errors) / offset

    sail_relative = worst(sail_points, sail_normals, sail, 0.01 * offset)
    fin_relative = worst(fin_points, fin_normals, fins, 0.002 * offset)
    # Two percent of the probe distance is the worst-case behaviour near the sail's
    # knife-edge tips; the median is three orders of magnitude better than that.
    assert sail_relative < 0.25
    assert fin_relative < 0.05


def test_the_normalised_gradient_is_unit_length() -> None:
    sail_points, sail_normals = _sail_surface_samples()
    fin_points, fin_normals = _fin_surface_samples()
    for field, points, normals in (
        (SuboffSailField(), sail_points, sail_normals),
        (SuboffFinsField(), fin_points, fin_normals),
    ):
        probes = [
            (p[0] + 5e-4 * n[0], p[1] + 5e-4 * n[1], p[2] + 5e-4 * n[2])
            for p, n in zip(points, normals, strict=True)
        ]
        numeric = field.gradient(probes)
        norms = [math.dist(gradient, (0.0, 0.0, 0.0)) for gradient in numeric]
        # The magnitude is the normalised gap, so |grad| <= 1 by construction and
        # approaches 1 where the surface is flat. It falls to 0.93 at the sail's
        # knife-edge tips, where the published section has an infinite slope.
        assert max(norms) <= 1.0 + 1e-9, field.kind
        assert sorted(norms)[len(norms) // 2] > 0.98, field.kind
        assert sum(norms) / len(norms) == pytest.approx(1.0, abs=0.02)


def test_the_flat_faces_report_the_exact_distance() -> None:
    """Where the surface is a plane the normalisation is exact, not first order."""
    feet = _feet()
    sail = SuboffSailField()
    width = _SAIL_HALF_WIDTH_FT / feet
    x = 3.4 / feet  # middle segment: a' = 0, so the side is a vertical plane
    for y_over in (1.2, 2.0, 5.0):
        point = (x, y_over * width, 0.4 * _SAIL_BODY_TOP_FT / feet)
        assert sail.signed_distance([point])[0] == pytest.approx((y_over - 1.0) * width, rel=1e-12)
    apex_x = 3.4 / feet
    top = (_SAIL_BODY_TOP_FT + 0.5 * _SAIL_HALF_WIDTH_FT) / feet
    assert sail.signed_distance([(apex_x, 0.0, top + 0.01)])[0] == pytest.approx(0.01, rel=1e-12)
    # ... and the fin's rim faces: the trailing edge is the plane x = _FIN_H, and the
    # tip is the cylinder r = r_outer, both exactly a distance away from outside.
    fins = SuboffFinsField()
    intercept, r_outer = _fin_chord_coefficients(1.0)
    radius = 0.4
    chord = -0.466308 * radius + intercept
    for step in (0.001, 0.01, 0.05):
        aft = (_FIN_TRAILING_FT + step * feet) / feet
        assert fins.signed_distance([(aft, 0.4 / feet, 0.0)])[0] == pytest.approx(step, rel=1e-12)
        assert fins.signed_distance([(aft, 0.0, 0.4 / feet)])[0] == pytest.approx(step, rel=1e-12)
    mid_chord = (_FIN_TRAILING_FT - 0.5 * chord) / feet
    for step in (0.001, 0.01):
        beyond = (r_outer + step * feet) / feet
        assert fins.signed_distance([(mid_chord, beyond, 0.0)])[0] == pytest.approx(step, rel=1e-12)
        assert fins.signed_distance([(mid_chord, 0.0, beyond)])[0] == pytest.approx(step, rel=1e-12)


# --------------------------------------------------------------------------- #
# Fields, edits, bounds
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    "factory",
    [SuboffSailField, SuboffFinsField, SuboffAppendagesField, SuboffFullField],
)
def test_invalid_parameters_and_queries_fail_loudly(factory: type) -> None:
    with pytest.raises(FieldError):
        factory(length=0.0)
    with pytest.raises(FieldError):
        factory(length=float("nan"))
    field = factory()
    with pytest.raises(FieldError):
        field.with_state(nonsense=1.0)
    with pytest.raises(FieldError):
        field.signed_distance([(0.0, 0.0)])
    with pytest.raises(FieldError):
        field.signed_distance([(0.0, 0.0, float("inf"))])
    assert field.signed_distance([]) == []
    assert field.contains([]) == []
    assert field.gradient([]) == []


def test_invalid_scales_fail_loudly() -> None:
    with pytest.raises(FieldError):
        SuboffSailField(scale=0.0)
    with pytest.raises(FieldError):
        SuboffFinsField(scale=-1.0)
    with pytest.raises(FieldError):
        SuboffAppendagesField(sail_scale=float("inf"))
    with pytest.raises(FieldError):
        SuboffFullField(fin_scale=0.0)


def test_bounds_are_the_published_extents() -> None:
    feet = _feet()
    sail = SuboffSailField()
    box = sail.bounds()
    assert box is not None
    assert box.minimum[0] == pytest.approx(_SAIL_X_START_FT / feet)
    assert box.maximum[0] == pytest.approx(_SAIL_X_END_FT / feet)
    assert box.maximum[1] == pytest.approx(_SAIL_HALF_WIDTH_FT / feet)
    assert box.maximum[2] == pytest.approx((_SAIL_BODY_TOP_FT + 0.5 * _SAIL_HALF_WIDTH_FT) / feet)
    fins = SuboffFinsField()
    box = fins.bounds()
    assert box is not None
    # The root chord is the longest, so the root leading edge is the forward bound.
    assert box.minimum[0] == pytest.approx((_FIN_TRAILING_FT - (-0.466308 * 0.075 + 0.88859)) / feet)
    assert box.minimum[1] == box.minimum[2] == pytest.approx(-0.825 / feet)
    assert box.maximum[1] == pytest.approx(0.825 / feet)
    full = SuboffFullField()
    box = full.bounds()
    assert box is not None
    assert box.minimum[0] == 0.0
    assert box.maximum[0] == pytest.approx(SUBOFF_LENGTH_M)
    assert box.maximum[2] == pytest.approx((_SAIL_BODY_TOP_FT + 0.5 * _SAIL_HALF_WIDTH_FT) / feet)


def test_the_bounds_contain_the_solid_they_claim() -> None:
    numpy = pytest.importorskip("numpy")
    generator = numpy.random.default_rng(11)
    field = SuboffFullField()
    box = field.bounds()
    assert box is not None
    count = 20_000
    points = numpy.column_stack(
        [
            generator.uniform(box.minimum[0] - 0.05, box.maximum[0] + 0.05, count),
            generator.uniform(box.minimum[1] - 0.05, box.maximum[1] + 0.05, count),
            generator.uniform(box.minimum[2] - 0.05, box.maximum[2] + 0.05, count),
        ]
    )
    inside = numpy.asarray(field.contains(points.tolist()))
    kept = points[inside]
    assert len(kept) > 100
    for axis in range(3):
        assert kept[:, axis].min() >= box.minimum[axis] - 1e-12
        assert kept[:, axis].max() <= box.maximum[axis] + 1e-12


def test_the_full_model_is_the_union_of_its_pieces() -> None:
    points = _to_metres([(x, y, z) for x, y, z, _, _ in TENSORLBM_POINTS_FT])
    full = SuboffFullField()
    hull = full.hull()
    appendages = full.appendages()
    union = [a or b for a, b in zip(hull.contains(points), appendages.contains(points), strict=True)]
    assert full.contains(points) == union
    values = full.signed_distance(points)
    hull_values = hull.signed_distance(points)
    appendage_values = appendages.signed_distance(points)
    assert values == [min(a, b) for a, b in zip(hull_values, appendage_values, strict=True)]
    metadata = full.metadata()
    assert metadata["kind"] == "suboff_full"
    assert metadata["quantity"] == IMPLICIT_FUNCTION
    assert metadata["units"] == "m"
    assert metadata["provenance"]["geometry_source"].startswith("DARPA SUBOFF")
    assert metadata["state"] == {"length": SUBOFF_LENGTH_M, "sail_scale": 1.0, "fin_scale": 1.0}
    assert full.parameters()["radius"] == pytest.approx(full.hull().max_radius)


def test_scaling_an_appendage_keeps_it_on_the_hull() -> None:
    """``scale`` grows the sail about its deck plane and the fins about their root."""
    feet = _feet()
    deck = _DECK_FT / feet
    sail = SuboffSailField()
    taller = sail.with_state(scale=2.0)
    # Height grows by `scale` about the deck plane, so the foot stays on the hull.
    assert taller.height == pytest.approx(deck + 2.0 * (sail.height - deck), rel=1e-12)
    assert sail.height == pytest.approx((_SAIL_BODY_TOP_FT + 0.5 * _SAIL_HALF_WIDTH_FT) / feet, rel=1e-12)
    assert taller.contains([(_SAIL_X_CENTER_FT / feet, 0.0, 0.99 * sail.height)]) == [True]
    assert taller.contains([(_SAIL_X_CENTER_FT / feet, 0.0, 0.99 * taller.height)]) == [True]
    assert sail.contains([(_SAIL_X_CENTER_FT / feet, 0.0, 0.99 * taller.height)]) == [False]
    # Width grows about the centreplane: the midbody half width doubles.
    assert float(_sail_half_width_array(3.4)) == pytest.approx(_SAIL_HALF_WIDTH_FT)
    sail_box = taller.bounds()
    assert sail_box is not None
    assert sail_box.maximum[1] == pytest.approx(2.0 * _SAIL_HALF_WIDTH_FT / feet, rel=1e-12)
    fins = SuboffFinsField()
    wider = fins.with_state(scale=1.25)
    assert wider.parameters()["tip_radius_m"] > fins.parameters()["tip_radius_m"]
    assert wider.parameters()["tip_radius_m"] == pytest.approx(
        (0.075 + 1.25 * (0.825 - 0.075)) / feet, rel=1e-12
    )
    tip_plane = (13.0 / feet, 0.8 / feet, 0.0)
    assert fins.contains([tip_plane]) == [True]
    assert wider.contains([tip_plane]) == [True]


def test_the_pure_python_path_matches_the_vectorised_one(monkeypatch: object) -> None:
    points = _to_metres(
        [
            (3.6, 0.0, 1.0),
            (3.6, 0.05, 1.53),
            (3.6, 0.2, 0.5),
            (12.9, 0.3, 0.05),
            (13.0, 0.0, 0.4),
            (3.0, 0.0, 1.0),
        ]
    )
    for factory in (SuboffSailField, SuboffFinsField, SuboffAppendagesField, SuboffFullField):
        field = factory()
        vectorised = field.signed_distance(points)
        monkeypatch.setattr(factory, "numpy_callable", lambda self: None)  # type: ignore[arg-type]
        assert field.numpy_callable() is None
        assert field.signed_distance(points) == vectorised


# --------------------------------------------------------------------------- #
# Export
# --------------------------------------------------------------------------- #


def _quadrature_sail_volume() -> float:
    """Sail volume from the published cross-section, independent of the field."""
    numpy = pytest.importorskip("numpy")
    feet = _feet()
    x = numpy.linspace(_SAIL_X_START_FT, _SAIL_X_END_FT, 400_001)
    width = _sail_half_width_array(x)
    area = 2.0 * width * _SAIL_BODY_TOP_FT + math.pi * width * width / 4.0
    return float(numpy.trapezoid(area, x)) / feet**3


def _quadrature_fin_volume() -> float:
    """Four fins' volume: ``2 h(s)`` integrated over each trapezoidal planform."""
    numpy = pytest.importorskip("numpy")
    feet = _feet()
    intercept, r_outer = _fin_chord_coefficients(1.0)
    s = numpy.linspace(0.0, 1.0, 400_001)
    integral = float(numpy.trapezoid(_naca_half_thickness_array(s), s))
    r = numpy.linspace(_FIN_R_INNER_FT, r_outer, 100_001)
    chord = -0.466308 * r + intercept
    return 4.0 * 2.0 * integral * float(numpy.trapezoid(chord, r)) / feet**3


def test_the_sail_exports_a_closed_mesh_at_the_right_volume() -> None:
    sampled = to_mesh(SuboffSailField(), (72, 12, 92), margin=2)
    assert sampled.report.watertight
    assert sampled.report.edge_manifold
    assert sampled.report.euler_characteristic == 2
    assert sampled.validation.errors == []
    reference = _quadrature_sail_volume()
    error = sampled.report.signed_volume / reference - 1.0
    assert abs(error) < 0.01, f"sail volume error {100 * error:+.3f} %"
    assert reference == pytest.approx(0.0091604303, abs=1e-8)


def test_the_fins_alone_are_a_ring_around_the_axis() -> None:
    """The blades start at r = 0.075 ft and are 0.1 ft thick, so they touch each other
    near the root and leave the axis column clear: one closed genus-1 shell."""
    sampled = to_mesh(SuboffFinsField(), (110, 64, 64), margin=2)
    assert sampled.report.watertight
    assert sampled.report.edge_manifold
    assert sampled.report.euler_characteristic == 0
    reference = _quadrature_fin_volume()
    error = sampled.report.signed_volume / reference - 1.0
    assert abs(error) < 0.02, f"fin volume error {100 * error:+.3f} %"
    assert reference == pytest.approx(0.0078571766, abs=1e-8)
    # The hole is real: the axis column is outside both blades ...
    assert SuboffFinsField().contains([(13.0 / _feet(), 0.02 / _feet(), 0.02 / _feet())]) == [False]


#: Sampling the full model is the expensive part of this file; the same field at the
#: same resolution is sampled once and asserted from several angles.
_MESH_CACHE: dict[tuple[int, ...], SampledField] = {}


def _full_mesh(resolution: tuple[int, int, int]) -> SampledField:
    if resolution not in _MESH_CACHE:
        _MESH_CACHE[resolution] = to_mesh(SuboffFullField(), resolution, margin=2)
    return _MESH_CACHE[resolution]


@pytest.mark.parametrize("resolution", [(160, 19, 28), (320, 38, 56)])
def test_the_full_model_exports_a_closed_solver_ready_mesh(resolution: tuple[int, int, int]) -> None:
    field = SuboffFullField()
    sampled = _full_mesh(resolution)
    assert sampled.report.watertight
    assert sampled.report.edge_manifold
    assert sampled.report.euler_characteristic == 2
    assert sampled.validation.errors == []
    volume_error = sampled.report.signed_volume / FULL_MODEL_VOLUME_M3 - 1.0
    assert abs(volume_error) < 0.01, f"volume error {100 * volume_error:+.3f} %"
    # The hull's own wetted area is a floor: the appendages add area, never remove it.
    hull_area = suboff_reference_geometry()["wetted_area"]
    assert sampled.report.surface_area > hull_area
    assert sampled.report.surface_area < 1.2 * hull_area
    # The mesh reaches the analytic bounds in z: the sail's cap apex.
    feet = _feet()
    apex = (_SAIL_BODY_TOP_FT + 0.5 * _SAIL_HALF_WIDTH_FT) / feet
    transverse_cell = (field.bounds().maximum[2] - field.bounds().minimum[2]) / resolution[2]  # type: ignore[union-attr]
    assert abs(sampled.report.bounds_max[2] - apex) < 2.0 * transverse_cell


def test_the_full_model_passes_a_handoff_record() -> None:
    record = _full_mesh((160, 19, 28)).handoff()
    assert record["units"] == "m"
    assert record["watertight"] is True
    assert record["euler_characteristic"] == 2
    assert record["validation"]["passed"] is True
    assert len(record["mesh_sha256"]) == 64
    assert record["signed_volume_m3"] > 0.69


def test_the_sampler_is_judged_against_the_formulas_not_against_itself() -> None:
    """Coarse and fine meshes bracket the quadrature volume from the same side."""
    coarse = _full_mesh((160, 19, 28))
    fine = _full_mesh((320, 38, 56))
    coarse_error = coarse.report.signed_volume / FULL_MODEL_VOLUME_M3 - 1.0
    fine_error = fine.report.signed_volume / FULL_MODEL_VOLUME_M3 - 1.0
    assert coarse_error < fine_error < 0.0
    assert abs(fine_error) < 0.005
    assert len(fine.mesh.faces) > len(coarse.mesh.faces)


def test_the_hull_profile_still_feeds_the_full_model() -> None:
    """The full model's hull half is the 0.5.0 field, unchanged."""
    hull = SuboffFullField().hull()
    reference = suboff_reference_geometry()
    for xi in (0.05, 0.5, 0.9):
        radius = float(suboff_radius_fraction_array(xi)) * hull.max_radius
        assert radius == pytest.approx(hull.radius_at(xi * hull.length), rel=1e-12)
    assert reference["volume"] == pytest.approx(0.699211, abs=1e-6)


def test_each_field_carries_its_own_record() -> None:
    """The pieces answer metadata on their own terms, not only through the union."""
    feet = _feet()
    sail = SuboffSailField()
    fins = SuboffFinsField()
    appendages = SuboffAppendagesField()
    for field in (sail, fins, appendages):
        metadata = field.metadata()
        assert metadata["quantity"] == IMPLICIT_FUNCTION
        assert metadata["units"] == "m"
        assert metadata["provenance"]["transcription"].startswith("TensorLBM suboff_cad")
        assert metadata["state"] == field.state()
        assert metadata["bounds"] == {
            key: list(value)
            for key, value in field.bounds().to_dict().items()  # type: ignore[union-attr]
        }
    assert sail.metadata()["kind"] == "suboff_sail"
    assert fins.metadata()["kind"] == "suboff_fins"
    assert sail.parameters()["half_width_m"] == pytest.approx(_SAIL_HALF_WIDTH_FT / feet, rel=1e-12)
    assert sail.parameters()["body_top_m"] == pytest.approx(_SAIL_BODY_TOP_FT / feet, rel=1e-12)
    assert fins.parameters()["root_radius_m"] == pytest.approx(_FIN_R_INNER_FT / feet, rel=1e-12)
    assert appendages.parameters()["pieces"] == [sail.parameters(), fins.parameters()]
    assert len(appendages.pieces()) == 2
    assert len(SuboffFullField().pieces()) == 2


def test_the_composites_delegate_and_edit_only_what_they_name() -> None:
    feet = _feet()
    hull_side = (2.0, 0.30, 0.0)  # outside the hull, far from any appendage
    sail_side = (_SAIL_X_CENTER_FT / feet, 0.0, 0.45)  # inside the sail
    points = [hull_side, sail_side]
    appendages = SuboffAppendagesField()
    full = SuboffFullField()
    # The gradient is the gradient of whichever piece is closer.
    sail_values = appendages.sail().signed_distance(points)
    fin_values = appendages.fins().signed_distance(points)
    gradients = appendages.gradient(points)
    sail_gradients = appendages.sail().gradient(points)
    fin_gradients = appendages.fins().gradient(points)
    for index, (s_value, f_value) in enumerate(zip(sail_values, fin_values, strict=True)):
        expected = sail_gradients[index] if s_value <= f_value else fin_gradients[index]
        assert gradients[index] == expected
    full_gradients = full.gradient(points)
    assert full_gradients[0] == pytest.approx(full.hull().gradient([hull_side])[0])
    assert full_gradients[1] == pytest.approx(full.appendages().gradient([sail_side])[0])
    # Edits are per-control and leave the others alone.
    bigger = appendages.with_state(sail_scale=2.0)
    assert bigger.state() == {"length": SUBOFF_LENGTH_M, "sail_scale": 2.0, "fin_scale": 1.0}
    assert bigger.fins().parameters() == appendages.fins().parameters()
    longer = full.with_state(length=8.0)
    assert longer.length == 8.0
    # The whole model scales: the hull keeps the published radius ratio, and the
    # appendages keep their published fractions of the length.
    assert longer.hull().max_radius / 8.0 == pytest.approx(full.hull().max_radius / full.length)
    assert (
        longer.appendages().sail().parameters()["half_width_m"]
        > (full.appendages().sail().parameters()["half_width_m"])
    )
    assert longer.appendages().length == 8.0
    assert full.with_state() == full


def test_the_scalar_path_walks_the_sail_entrance_segment(monkeypatch: object) -> None:
    """The pure-Python path has its own slope code; hit the entrance branch of it."""
    feet = _feet()
    points = [
        (3.05 / feet, 0.0, 1.0 / feet),
        (3.05 / feet, 0.05 / feet, 1.0 / feet),
        (4.20 / feet, 0.0, 1.0 / feet),
        (3.40 / feet, 0.10 / feet, 1.55 / feet),
    ]
    sail = SuboffSailField()
    vectorised = sail.signed_distance(points)
    monkeypatch.setattr(SuboffSailField, "numpy_callable", lambda self: None)  # type: ignore[arg-type]
    assert sail.numpy_callable() is None
    assert sail.signed_distance(points) == vectorised
    # The entrance slope is finite at 3.05 ft and unbounded only at the leading tip.
    assert math.isfinite(float(_sail_half_width_slope_array(3.05)))
    assert float(_sail_half_width_slope_array(3.05)) > 0.0
