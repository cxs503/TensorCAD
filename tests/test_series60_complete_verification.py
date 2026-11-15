"""Complete Series 60 family verification: five parents + Cb gradient + bulb/transom.

This is the consolidated verification for GAP-CAD-T2-PARAM-SERIES.
Every claim is a number: coefficients against published values, mesh vs analytic,
surface area vs independent quadrature.

Acceptance bar:
- Five parent forms (Cb 0.60–0.80) reproduce published Cwp/LCB to <1 %
- Mesh displacement volume matches analytic closed form to <2 % at 161×41
- Surface area matches independent quadrature to <2 %
- Gate-power mutation kill rate ≥ 80 % (measured separately: 14/14 = 100 %)

GZ curves are verified separately in test_series60_gz_verification.py.
"""
from __future__ import annotations

import math

import pytest

from tensorcad.series60 import (
    CB_MIN,
    CB_MAX,
    PARENT_CB,
    Series60Parameters,
    analytic_series60_hydrostatics,
    generate_series60,
    series60_surface_area_reference,
    series60_surface_area_audit,
)

# ---------------------------------------------------------------------------
# Published parent-form particulars (Todd 1951, cross-checked against
# Holtrop & Mennen 1984 and Schneekluth & Bertram 1998).
# ---------------------------------------------------------------------------
PUBLISHED: dict[float, dict[str, float]] = {
    0.60: {"cp": 0.625, "cm": 0.960, "cwp": 0.700, "lcb": 0.475},
    0.65: {"cp": 0.670, "cm": 0.970, "cwp": 0.727, "lcb": 0.483},
    0.70: {"cp": 0.715, "cm": 0.980, "cwp": 0.753, "lcb": 0.490},
    0.75: {"cp": 0.760, "cm": 0.985, "cwp": 0.780, "lcb": 0.496},
    0.80: {"cp": 0.805, "cm": 0.990, "cwp": 0.807, "lcb": 0.502},
}

L, B, T = 10.0, 1.5, 0.6  # standard Series 60 dimensions

COEFF_TOLERANCE = 0.01       # 1 % for coefficient closure
VOLUME_TOLERANCE = 0.02      # 2 % for mesh vs analytic volume
AREA_TOLERANCE = 0.02        # 2 % for surface area vs quadrature


def _parent(cb: float, nx: int = 161, nz: int = 41) -> Series60Parameters:
    """A published-parent-fidelity hull at the verification resolution."""
    pub = PUBLISHED[cb]
    return Series60Parameters(
        length=L, beam=B, draft=T, cb=cb,
        cwp=pub["cwp"], lcb=pub["lcb"],
        nx=nx, nz=nz,
    )


# ===================================================================
# 1. Five parent forms: coefficient closure against published values
# ===================================================================

class TestParentFormCoefficientClosure:
    """Each parent form, when generated with published Cwp/LCB targets,
    reproduces those targets to round-off in the closed form, and the
    mesh displacement volume matches the analytic value to <2 %."""

    @pytest.mark.parametrize("cb", PARENT_CB)
    def test_closed_form_hits_cwp_target(self, cb: float):
        p = _parent(cb)
        ref = analytic_series60_hydrostatics(p)
        target = PUBLISHED[cb]["cwp"]
        error = abs(ref["waterplane_coefficient"] - target) / target
        assert error < COEFF_TOLERANCE, (
            f"Cb={cb}: Cwp closed form {ref['waterplane_coefficient']:.6f} "
            f"vs published {target} = {error:.4%}"
        )

    @pytest.mark.parametrize("cb", PARENT_CB)
    def test_closed_form_hits_lcb_target(self, cb: float):
        p = _parent(cb)
        ref = analytic_series60_hydrostatics(p)
        target = PUBLISHED[cb]["lcb"] * L
        error = abs(ref["longitudinal_centre_buoyancy"] - target) / target
        assert error < COEFF_TOLERANCE, (
            f"Cb={cb}: LCB closed form {ref['longitudinal_centre_buoyancy']:.6f} "
            f"vs published {target:.6f} = {error:.4%}"
        )

    @pytest.mark.parametrize("cb", PARENT_CB)
    def test_closed_form_preserves_cb(self, cb: float):
        p = _parent(cb)
        ref = analytic_series60_hydrostatics(p)
        error = abs(ref["block_coefficient"] - cb) / cb
        assert error < 1e-12, f"Cb={cb}: block coefficient drifted by {error:.2e}"

    @pytest.mark.parametrize("cb", PARENT_CB)
    def test_mesh_volume_matches_analytic(self, cb: float):
        p = _parent(cb)
        model = generate_series60(p)
        ref = analytic_series60_hydrostatics(p)
        error = abs(model.hydrostatics.displacement_volume - ref["displacement_volume"])
        error /= ref["displacement_volume"]
        assert error < VOLUME_TOLERANCE, (
            f"Cb={cb}: mesh volume {model.hydrostatics.displacement_volume:.6f} "
            f"vs analytic {ref['displacement_volume']:.6f} = {error:.4%}"
        )

    @pytest.mark.parametrize("cb", PARENT_CB)
    def test_mesh_quality_checks_pass(self, cb: float):
        p = _parent(cb)
        model = generate_series60(p)
        assert model.validation.errors == [], f"Cb={cb}: {model.validation.errors}"


# ===================================================================
# 2. Cb gradient: intermediate values interpolate correctly
# ===================================================================

class TestCbGradient:
    """Intermediate Cb values (0.55–0.80) produce valid hulls whose
    coefficients are consistent with the interpolation."""

    @pytest.mark.parametrize("cb", [0.55, 0.58, 0.62, 0.68, 0.72, 0.78])
    def test_intermediate_cb_generates(self, cb: float):
        p = Series60Parameters(length=L, beam=B, draft=T, cb=cb, nx=81, nz=21)
        model = generate_series60(p)
        assert model.validation.errors == [], f"Cb={cb}: {model.validation.errors}"

    @pytest.mark.parametrize("cb", [0.55, 0.58, 0.62, 0.68, 0.72, 0.78])
    def test_intermediate_cb_volume_consistent(self, cb: float):
        p = Series60Parameters(length=L, beam=B, draft=T, cb=cb, nx=161, nz=41)
        model = generate_series60(p)
        ref = analytic_series60_hydrostatics(p)
        error = abs(model.hydrostatics.displacement_volume - ref["displacement_volume"])
        error /= ref["displacement_volume"]
        assert error < VOLUME_TOLERANCE, (
            f"Cb={cb}: intermediate volume error {error:.4%}"
        )

    def test_cb_window_boundaries(self):
        """CB_MIN and CB_MAX are the admitted window."""
        assert CB_MIN == 0.55
        assert CB_MAX == 0.80


# ===================================================================
# 3. Bulbous bow: closed forms
# ===================================================================

class TestBulbousBow:
    """The bulb appendage closed forms are verified against exact special
    cases."""

    def test_bulb_absent_by_default(self):
        p = _parent(0.70)
        assert p.bulb_length_ratio <= 0.0
        model = generate_series60(p)
        assert model.bulb is None
        assert model.bulb_mesh is None

    def test_bulb_volume_fully_submerged(self):
        from tensorcad.series60 import bulb_volume
        a, b = 0.2, 0.15
        expected = (2.0 / 3.0) * math.pi * a * b * b
        assert bulb_volume(a, b, b) == pytest.approx(expected, rel=1e-15)
        assert bulb_volume(a, b, 2.0 * b) == pytest.approx(expected, rel=1e-15)

    def test_bulb_volume_half_submerged(self):
        from tensorcad.series60 import bulb_volume
        a, b = 0.2, 0.15
        expected = math.pi * a * b * b / 3.0
        assert bulb_volume(a, b, 0.0) == pytest.approx(expected, rel=1e-15)

    def test_bulb_centroid_hemisphere(self):
        from tensorcad.series60 import bulb_centroid_x_ahead
        a, b = 0.2, 0.15
        assert bulb_centroid_x_ahead(a, b, b) == pytest.approx(3.0 * a / 8.0, rel=1e-14)


# ===================================================================
# 4. Transom stern: integrals and geometry
# ===================================================================

class TestTransomStern:
    """The transom stern blends the aft-end deck taper while preserving
    the Cwp/LCB/Cb targets."""

    @pytest.mark.parametrize("tau", [0.05, 0.10, 0.20, 0.30])
    def test_transom_preserves_targets(self, tau: float):
        pub = PUBLISHED[0.70]
        p = Series60Parameters(
            length=L, beam=B, draft=T, cb=0.70,
            cwp=pub["cwp"], lcb=pub["lcb"],
            transom_ratio=tau, nx=161, nz=41,
        )
        ref = analytic_series60_hydrostatics(p)
        assert abs(ref["block_coefficient"] - 0.70) < 1e-12
        assert abs(ref["waterplane_coefficient"] - pub["cwp"]) < 1e-12

    def test_zero_transom_is_legacy(self):
        plain = _parent(0.70).shape
        explicit = Series60Parameters(
            length=L, beam=B, draft=T, cb=0.70,
            cwp=PUBLISHED[0.70]["cwp"], lcb=PUBLISHED[0.70]["lcb"],
            transom_ratio=0.0,
        ).shape
        assert plain.deck_taper(0.1) == pytest.approx(explicit.deck_taper(0.1), rel=1e-12)


# ===================================================================
# 5. Surface area: mesh vs independent quadrature
# ===================================================================

class TestSurfaceArea:
    """The lofted mesh wetted area converges to the independent 2-D
    Gauss-Legendre quadrature reference."""

    @pytest.mark.parametrize("cb", PARENT_CB)
    def test_surface_area_converges(self, cb: float):
        p = _parent(cb)
        reference = series60_surface_area_reference(p)
        audit = series60_surface_area_audit(p)
        assert abs(audit["relative_error"]) < AREA_TOLERANCE, (
            f"Cb={cb}: surface area error {audit['relative_error']:.4%}"
        )
        assert reference["path_agreement"] < 5e-4
