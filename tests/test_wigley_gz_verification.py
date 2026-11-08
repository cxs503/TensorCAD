"""Wigley hull GZ curve verification against analytical GM.

The Wigley parabolic hull is the standard benchmark for wave resistance and
seakeeping. This module verifies that the GZ integrator produces a righting-arm
curve whose small-angle slope matches the metacentric height GM from waterplane
moments -- two independent integrations of the same hull.

The Wigley hull has no freeboard by default (deck at waterline), so it must be
extended with ``with_freeboard`` before computing GZ. The waterplane moments
are computed on the original mesh (the freeboard band is vertical and does not
contribute to the waterplane).
"""
import math
import pytest

from tensorcad.wigley import WigleyParameters, generate_wigley
from tensorcad.gz import with_freeboard
from tensorcad.righting_arm import righting_arm_curve, heeled_displacement
from tensorcad.stability import waterplane_moments


FREEBOARD = 0.2
KG = 0.3
WATERLINE = 0.0


def wigley_model(nx: int = 81, nz: int = 21):
    """Generate a Wigley hull model."""
    p = WigleyParameters(length=10.0, beam=1.0, draft=0.625, nx=nx, nz=nz)
    return generate_wigley(p)


def compute_gm_from_original_mesh(model):
    """Compute GM from the original mesh (before freeboard).
    
    The waterplane is the same with or without freeboard, since the freeboard
    band is vertical above z=0.
    """
    mesh = model.mesh
    state = heeled_displacement(mesh, heel_deg=0.0, waterline_z=WATERLINE)
    section = waterplane_moments(mesh, WATERLINE)
    keel_z = mesh.bounding_box().minimum[2]
    kb = state.centre_of_buoyancy[2] - keel_z
    bm = section.transverse_second_moment / state.volume
    return kb + bm - KG


class TestWigleyGZSmallAngleSlope:
    """The slope of GZ at φ→0 must equal GM from waterplane moments."""

    def test_slope_matches_gm_at_5_degrees(self):
        """GZ/sin(φ) at 5° should match GM to within 0.5%."""
        model = wigley_model()
        gm = compute_gm_from_original_mesh(model)
        
        fb_mesh = with_freeboard(model.mesh, freeboard=FREEBOARD, waterline_z=WATERLINE)
        curve = righting_arm_curve(
            fb_mesh, kg=KG, waterline_z=WATERLINE,
            heels_deg=(0.0, 5.0)
        )
        
        phi_5 = math.radians(5.0)
        gz_5 = curve.points[1].righting_arm
        slope = gz_5 / math.sin(phi_5)
        
        assert slope == pytest.approx(gm, rel=0.005), (
            f"GZ slope {slope:.6f} != GM {gm:.6f} (ratio {slope/gm:.4f})"
        )

    def test_slope_matches_gm_at_10_degrees(self):
        """GZ/sin(φ) at 10° should match GM to within 1.5%."""
        model = wigley_model()
        gm = compute_gm_from_original_mesh(model)
        
        fb_mesh = with_freeboard(model.mesh, freeboard=FREEBOARD, waterline_z=WATERLINE)
        curve = righting_arm_curve(
            fb_mesh, kg=KG, waterline_z=WATERLINE,
            heels_deg=(0.0, 10.0)
        )
        
        phi_10 = math.radians(10.0)
        gz_10 = curve.points[1].righting_arm
        slope = gz_10 / math.sin(phi_10)
        
        # At 10°, the wall-sided tan²φ term contributes ~1.5%
        assert slope == pytest.approx(gm, rel=0.015), (
            f"GZ slope {slope:.6f} != GM {gm:.6f} (ratio {slope/gm:.4f})"
        )


class TestWigleyVolumeConservation:
    """Displaced volume must be constant at each heel angle."""

    def test_volume_conserved_to_machine_precision(self):
        """Volume error must be < 1e-10 relative."""
        model = wigley_model()
        fb_mesh = with_freeboard(model.mesh, freeboard=FREEBOARD, waterline_z=WATERLINE)
        
        curve = righting_arm_curve(
            fb_mesh, kg=KG, waterline_z=WATERLINE,
            heels_deg=tuple(range(0, 46, 5))
        )
        
        V_upright = curve.points[0].displacement_volume
        
        for pt in curve.points:
            rel_error = abs(pt.displacement_volume - V_upright) / V_upright
            assert rel_error < 1e-10, (
                f"Volume not conserved at heel={pt.heel_deg}°: "
                f"rel error {rel_error:.2e}"
            )


class TestWigleyGZPhysicalProperties:
    """GZ curve must have expected physical properties."""

    def test_gz_vanishes_upright(self):
        """GZ(0) = 0."""
        model = wigley_model()
        fb_mesh = with_freeboard(model.mesh, freeboard=FREEBOARD, waterline_z=WATERLINE)
        curve = righting_arm_curve(
            fb_mesh, kg=KG, waterline_z=WATERLINE,
            heels_deg=(0.0, 10.0)
        )
        assert curve.points[0].righting_arm == pytest.approx(0.0, abs=1e-12)

    def test_gz_positive_for_small_heel(self):
        """GZ > 0 for small heel angles (stable equilibrium)."""
        model = wigley_model()
        fb_mesh = with_freeboard(model.mesh, freeboard=FREEBOARD, waterline_z=WATERLINE)
        curve = righting_arm_curve(
            fb_mesh, kg=KG, waterline_z=WATERLINE,
            heels_deg=tuple(range(0, 31, 5))
        )
        for pt in curve.points:
            if pt.heel_deg > 0:
                assert pt.righting_arm > 0.0, (
                    f"GZ negative at heel={pt.heel_deg}°: {pt.righting_arm:.6f}"
                )

    def test_gz_increases_monotonically_initially(self):
        """GZ should increase for small heel angles."""
        model = wigley_model()
        fb_mesh = with_freeboard(model.mesh, freeboard=FREEBOARD, waterline_z=WATERLINE)
        curve = righting_arm_curve(
            fb_mesh, kg=KG, waterline_z=WATERLINE,
            heels_deg=tuple(range(0, 31, 5))
        )
        for i in range(1, len(curve.points)):
            assert curve.points[i].righting_arm > curve.points[i-1].righting_arm, (
                f"GZ not increasing at heel={curve.points[i].heel_deg}°"
            )


class TestWigleyAnalyticHydrostatics:
    """Wigley hull hydrostatics match analytical closed forms."""

    def test_volume_matches_analytic(self):
        """Displaced volume matches Cb·L·B·T = 4LBT/9."""
        model = wigley_model(nx=161, nz=41)
        L, B, T = 10.0, 1.0, 0.625
        V_analytic = 4.0 * L * B * T / 9.0
        V_mesh = model.hydrostatics.displacement_volume
        rel_error = abs(V_mesh - V_analytic) / V_analytic
        assert rel_error < 0.002, (
            f"Volume error {rel_error:.4f}: mesh {V_mesh:.6f} vs analytic {V_analytic:.6f}"
        )

    def test_waterplane_area_matches_analytic(self):
        """Waterplane area matches 2LB/3."""
        model = wigley_model(nx=161, nz=41)
        L, B = 10.0, 1.0
        Awp_analytic = 2.0 * L * B / 3.0
        Awp_mesh = model.hydrostatics.waterplane_area
        rel_error = abs(Awp_mesh - Awp_analytic) / Awp_analytic
        assert rel_error < 0.005, (
            f"Waterplane area error {rel_error:.4f}: mesh {Awp_mesh:.6f} vs analytic {Awp_analytic:.6f}"
        )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
