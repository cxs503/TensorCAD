"""Tests for :mod:`tensorcad.box_inclined` -- exact barge hydrostatics at heel + trim.

Two independent routes are used as references, and neither is the module's own algorithm:

* the wall-sided closed form ``GZ = sin(phi)(GM + BM/2 tan^2 phi)`` (upright window only),
* a **2-D section oracle** written here: the rectangular section is clipped by a horizontal
  line and its area/centroid are integrated from the resulting polygon, i.e. a dimension
  below the module's 3-D polyhedron clip.  The two routes share no code.
"""
from __future__ import annotations

import math

import pytest

from tensorcad.box_inclined import (
    BoxInclinedError,
    box_inclined_curve,
    box_inclined_gz,
    box_inclined_hydrostatics,
    submerged_box,
    waterline_normal,
)

LENGTH, BEAM, DEPTH, DRAFT = 10.0, 2.0, 1.2, 0.8
VOLUME = LENGTH * BEAM * DRAFT
KG = 0.55


# ----------------------------------------------------------------------------------
# 2-D section oracle (independent route: one dimension lower, different geometry code)
# ----------------------------------------------------------------------------------
def section_oracle(heel_deg: float, *, beam=BEAM, draft=DRAFT, depth=DEPTH):
    """Submerged 2-D section in **body** coordinates: ``(area, y_centroid, z_centroid)``.

    The rectangle is clipped directly by the body-frame waterline
    ``-sin(phi) y + cos(phi) z <= c``; the offset ``c`` is bisected until the submerged area
    equals the prismatic displacement per unit length.  One dimension below the module's
    3-D polyhedron, and no shared code -- though it does share the half-space idea, so the
    agreement it demonstrates is "two implementations of the same cut", not "two methods".
    """
    phi = math.radians(heel_deg)
    ny, nz = -math.sin(phi), math.cos(phi)
    area2d = beam * draft
    polygon = [(-beam / 2.0, 0.0), (beam / 2.0, 0.0), (beam / 2.0, depth), (-beam / 2.0, depth)]

    def clipped(c):
        below = []
        count = len(polygon)
        for index in range(count):
            y0, z0 = polygon[index]
            y1, z1 = polygon[(index + 1) % count]
            d0 = ny * y0 + nz * z0 - c
            d1 = ny * y1 + nz * z1 - c
            if d0 <= 0.0:
                below.append((y0, z0))
            if (d0 < 0.0 < d1) or (d1 < 0.0 < d0):
                t = d0 / (d0 - d1)
                below.append((y0 + t * (y1 - y0), z0 + t * (z1 - z0)))
        return below

    def props(c):
        poly = clipped(c)
        if len(poly) < 3:
            return 0.0, 0.0, 0.0
        a = cy = cz = 0.0
        for index in range(len(poly)):
            y0, z0 = poly[index]
            y1, z1 = poly[(index + 1) % len(poly)]
            cross = y0 * z1 - y1 * z0
            a += cross
            cy += (y0 + y1) * cross
            cz += (z0 + z1) * cross
        a *= 0.5
        return abs(a), cy / (6.0 * a), cz / (6.0 * a)

    lo, hi = -2.0 * (beam + depth), 2.0 * (beam + depth)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if props(mid)[0] < area2d:
            lo = mid
        else:
            hi = mid
    area, cy, cz = props(0.5 * (lo + hi))
    return area, cy, cz


def wall_sided_gz(heel_deg, *, beam=BEAM, draft=DRAFT, kg=KG):
    bm = beam * beam / (12.0 * draft)
    gm = 0.5 * draft + bm - kg
    phi = math.radians(heel_deg)
    return math.sin(phi) * (gm + bm / 2.0 * math.tan(phi) ** 2)


# ----------------------------------------------------------------------------------
# upright limits
# ----------------------------------------------------------------------------------
class TestUprightLimits:
    def test_volume_and_buoyancy_centre(self):
        hydro = box_inclined_hydrostatics(LENGTH, BEAM, DEPTH, volume_m3=VOLUME)
        assert hydro.volume_m3 == pytest.approx(VOLUME, rel=1e-14)
        assert hydro.lcb == pytest.approx(0.0, abs=1e-12)
        assert hydro.tcb == pytest.approx(0.0, abs=1e-12)
        assert hydro.vcb == pytest.approx(DRAFT / 2.0, rel=1e-14)

    def test_waterplane_moments_match_the_rectangle(self):
        hydro = box_inclined_hydrostatics(LENGTH, BEAM, DEPTH, volume_m3=VOLUME)
        assert hydro.waterplane_area_m2 == pytest.approx(LENGTH * BEAM, rel=1e-14)
        assert hydro.i_transverse_m4 == pytest.approx(LENGTH * BEAM**3 / 12.0, rel=1e-14)
        assert hydro.i_longitudinal_m4 == pytest.approx(BEAM * LENGTH**3 / 12.0, rel=1e-13)
        assert hydro.bm_transverse == pytest.approx(BEAM**2 / (12.0 * DRAFT), rel=1e-14)
        assert hydro.bm_longitudinal == pytest.approx(LENGTH**2 / (12.0 * DRAFT), rel=1e-14)

    def test_wetted_area_is_bottom_plus_four_sides(self):
        hydro = box_inclined_hydrostatics(LENGTH, BEAM, DEPTH, volume_m3=VOLUME)
        expected = LENGTH * BEAM + 2.0 * LENGTH * DRAFT + 2.0 * BEAM * DRAFT
        assert hydro.wetted_area_m2 == pytest.approx(expected, rel=1e-13)


# ----------------------------------------------------------------------------------
# trim only: closed forms derived from the plane cut, independent of the polyhedron route
# ----------------------------------------------------------------------------------
class TestTrimOnly:
    @pytest.mark.parametrize("trim_deg", [1.0, 2.0, 3.0, 4.0])
    def test_lcb_matches_the_analytic_trim_moment(self, trim_deg):
        # submerged thickness a + x tan(tau) with a = V/(L B): LCB = L^2 tan(tau) / (12 T)
        hydro = box_inclined_hydrostatics(LENGTH, BEAM, DEPTH, volume_m3=VOLUME, trim_deg=trim_deg)
        expected = LENGTH**2 * math.tan(math.radians(trim_deg)) / (12.0 * DRAFT)
        assert hydro.volume_m3 == pytest.approx(VOLUME, rel=1e-14)
        assert hydro.lcb == pytest.approx(expected, rel=1e-9)
        assert hydro.tcb == pytest.approx(0.0, abs=1e-12)

    def test_waterplane_stretches_by_one_over_cos_trim(self):
        trim_deg = 4.0
        hydro = box_inclined_hydrostatics(LENGTH, BEAM, DEPTH, volume_m3=VOLUME, trim_deg=trim_deg)
        assert hydro.waterplane_area_m2 == pytest.approx(
            LENGTH * BEAM / math.cos(math.radians(trim_deg)), rel=1e-13
        )
        assert hydro.i_transverse_m4 == pytest.approx(
            LENGTH / math.cos(math.radians(trim_deg)) * BEAM**3 / 12.0, rel=1e-13
        )

    def test_trim_shifts_no_transverse_centre(self):
        hydro = box_inclined_hydrostatics(LENGTH, BEAM, DEPTH, volume_m3=VOLUME, trim_deg=5.0)
        assert hydro.tcb == pytest.approx(0.0, abs=1e-12)


# ----------------------------------------------------------------------------------
# heel: wall-sided window and the independent 2-D section route
# ----------------------------------------------------------------------------------
class TestHeel:
    @pytest.mark.parametrize("heel_deg", [1.0, 2.0, 5.0, 10.0, 15.0])
    def test_inside_the_wall_sided_window(self, heel_deg):
        gz = box_inclined_gz(
            LENGTH, BEAM, DEPTH, kg=KG, volume_m3=VOLUME, heel_deg=heel_deg
        )
        assert gz == pytest.approx(wall_sided_gz(heel_deg), rel=1e-12, abs=1e-15)

    @pytest.mark.parametrize("heel_deg", [10.0, 20.0, 30.0, 40.0, 55.0, 70.0])
    def test_matches_the_2d_section_oracle(self, heel_deg):
        hydro = box_inclined_hydrostatics(
            LENGTH, BEAM, DEPTH, volume_m3=VOLUME, heel_deg=heel_deg
        )
        area, y2d, z2d = section_oracle(heel_deg)
        assert hydro.volume_m3 == pytest.approx(VOLUME, rel=1e-13)
        # the prismatic section centroid is the 3-D buoyancy centre's (y, z)
        assert hydro.tcb == pytest.approx(y2d, rel=0, abs=2e-9)
        assert hydro.vcb == pytest.approx(z2d, rel=0, abs=2e-9)
        assert area == pytest.approx(BEAM * DRAFT, rel=1e-12)

    def test_the_wall_sided_formula_is_not_vacuous_here(self):
        """Past its window the closed form departs from the exact arm, so the check bites."""
        heel = 40.0
        exact = box_inclined_gz(LENGTH, BEAM, DEPTH, kg=KG, volume_m3=VOLUME, heel_deg=heel)
        assert abs(exact - wall_sided_gz(heel)) > 0.02

    def test_mirrored_heel_mirrors_the_buoyancy_centre(self):
        starboard = box_inclined_hydrostatics(
            LENGTH, BEAM, DEPTH, volume_m3=VOLUME, heel_deg=25.0
        )
        port = box_inclined_hydrostatics(
            LENGTH, BEAM, DEPTH, volume_m3=VOLUME, heel_deg=-25.0
        )
        assert starboard.tcb == pytest.approx(-port.tcb, abs=1e-12)
        assert starboard.vcb == pytest.approx(port.vcb, rel=1e-12)

    def test_heel_and_trim_compose_without_losing_the_volume(self):
        hydro = box_inclined_hydrostatics(
            LENGTH, BEAM, DEPTH, volume_m3=VOLUME, heel_deg=18.0, trim_deg=3.0
        )
        assert hydro.volume_m3 == pytest.approx(VOLUME, rel=1e-13)
        assert hydro.lcb > 0.0  # bow down
        assert hydro.tcb > 0.0  # starboard down
        # both effects are first order and independent at this size
        trim_only = box_inclined_hydrostatics(
            LENGTH, BEAM, DEPTH, volume_m3=VOLUME, trim_deg=3.0
        )
        # 18 deg of heel re-cuts the waterplane, so the trim-only LCB is a 1.4 % proxy
        assert hydro.lcb == pytest.approx(trim_only.lcb, rel=2e-2)


# ----------------------------------------------------------------------------------
# the arm: sign, peak and monotone growth at small heel
# ----------------------------------------------------------------------------------
class TestRightingArm:
    def test_upright_arm_is_zero_and_slope_is_gm(self):
        bm = BEAM**2 / (12.0 * DRAFT)
        gm = 0.5 * DRAFT + bm - KG
        small = box_inclined_gz(LENGTH, BEAM, DEPTH, kg=KG, volume_m3=VOLUME, heel_deg=1e-4)
        assert box_inclined_gz(LENGTH, BEAM, DEPTH, kg=KG, volume_m3=VOLUME, heel_deg=0.0) == pytest.approx(0.0, abs=1e-12)
        assert small / math.radians(1e-4) == pytest.approx(gm, rel=1e-6)

    def test_curve_rises_then_falls(self):
        angles = tuple(range(0, 91, 5))
        curve = box_inclined_curve(
            LENGTH, BEAM, DEPTH, kg=KG, volume_m3=VOLUME, angles_deg=angles
        )
        peak = max(range(len(curve)), key=lambda i: curve[i])
        assert 0 < peak < len(curve) - 1
        assert curve[peak] > curve[peak - 1] and curve[peak] >= curve[peak + 1]

    def test_trimming_moves_the_arm(self):
        untrimmed = box_inclined_gz(LENGTH, BEAM, DEPTH, kg=KG, volume_m3=VOLUME, heel_deg=20.0)
        trimmed = box_inclined_gz(
            LENGTH, BEAM, DEPTH, kg=KG, volume_m3=VOLUME, heel_deg=20.0, trim_deg=4.0
        )
        assert trimmed != pytest.approx(untrimmed, rel=1e-6)


# ----------------------------------------------------------------------------------
# the polyhedron itself: closure and refusals
# ----------------------------------------------------------------------------------
class TestPolyhedron:
    def test_clip_is_a_closed_convex_solid(self):
        hydro = submerged_box(
            LENGTH, BEAM, DEPTH, heel_deg=20.0, trim_deg=3.0, waterline_offset=0.4
        )
        assert hydro.polyhedron_volume_m3 == pytest.approx(hydro.volume_m3, rel=1e-14)
        assert hydro.faces >= 4
        # every vertex coordinate stays inside the box
        for point in hydro.waterplane_polygon:
            assert -LENGTH / 2 - 1e-9 <= point[0] <= LENGTH / 2 + 1e-9
            assert -BEAM / 2 - 1e-9 <= point[1] <= BEAM / 2 + 1e-9
            assert -1e-9 <= point[2] <= DEPTH + 1e-9

    def test_normal_limits(self):
        assert waterline_normal(0.0, 0.0) == pytest.approx((0.0, 0.0, 1.0), abs=1e-15)
        phi = math.radians(20.0)
        assert waterline_normal(20.0, 0.0) == pytest.approx((0.0, -math.sin(phi), math.cos(phi)), rel=1e-15)
        tau = math.radians(3.0)
        assert waterline_normal(0.0, 3.0) == pytest.approx((-math.sin(tau), 0.0, math.cos(tau)), rel=1e-15)

    def test_impossible_requests_are_refused(self):
        with pytest.raises(BoxInclinedError):
            box_inclined_hydrostatics(LENGTH, BEAM, DEPTH, volume_m3=LENGTH * BEAM * DEPTH * 1.5)
        with pytest.raises(BoxInclinedError):
            box_inclined_hydrostatics(LENGTH, BEAM, DEPTH, volume_m3=-1.0)
        with pytest.raises(BoxInclinedError):
            box_inclined_hydrostatics(0.0, BEAM, DEPTH, volume_m3=VOLUME)
        with pytest.raises(BoxInclinedError):
            submerged_box(LENGTH, BEAM, DEPTH, waterline_offset=-5.0)
        with pytest.raises(BoxInclinedError):
            box_inclined_hydrostatics(LENGTH, BEAM, DEPTH, volume_m3=VOLUME, heel_deg=float("nan"))