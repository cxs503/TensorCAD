"""Series 60 GZ verification: righting-arm curves for all five parent forms.

Verifies physical properties (GZ(0)=0, volume conservation, positive stability),
small-angle theory (dGZ/dφ|₀ = GM), and wall-sided formula agreement.

The Series 60 hull, like the Wigley hull, has no freeboard by default (deck at
waterline). The ``with_freeboard`` function must be called before computing GZ
curves, otherwise the waterline leaves the modelled surface on the dipping side
and the small-angle slope does not match GM.
"""
import math
import pytest
from tensorcad.series60 import Series60Parameters, generate_series60
from tensorcad.stability import initial_stability
from tensorcad.righting_arm import righting_arm_curve
from tensorcad.gz import with_freeboard


CB_VALUES = [0.60, 0.65, 0.70, 0.75, 0.80]
L, B, T = 10.0, 1.5, 0.6
KG = T * 0.4
NX, NZ = 41, 13  # coarser grid for test speed
FREEBOARD = 0.2
HEELS_DEG = tuple(range(0, 31, 5))


def wall_sided_gz(phi_rad: float, gm: float, bm: float) -> float:
    """Wall-sided formula: GZ = sin(phi)(GM + BM/2 tan^2(phi))."""
    return math.sin(phi_rad) * (gm + 0.5 * bm * math.tan(phi_rad) ** 2)


@pytest.mark.parametrize("cb", CB_VALUES)
def test_physical_properties(cb):
    """GZ(0)=0, volume conservation, GZ>0 for small heel."""
    p = Series60Parameters(length=L, beam=B, draft=T, cb=cb, nx=NX, nz=NZ)
    model = generate_series60(p)
    fb_mesh = with_freeboard(model.mesh, freeboard=FREEBOARD, waterline_z=0.0)
    curve = righting_arm_curve(fb_mesh, kg=KG, heels_deg=HEELS_DEG, waterline_z=0.0)
    
    # GZ(0) = 0
    assert abs(curve.points[0].righting_arm) < 1e-10
    
    # Volume conservation
    V0 = curve.points[0].displacement_volume
    for pt in curve.points:
        assert abs(pt.displacement_volume - V0) / V0 < 1e-8
    
    # GZ > 0 for small heel (stable)
    for pt in curve.points[1:4]:
        assert pt.righting_arm > 0


def test_small_angle_slope_cb070():
    """dGZ/dphi|_{phi=0} = GM for Cb=0.70."""
    p = Series60Parameters(length=L, beam=B, draft=T, cb=0.70, nx=NX, nz=NZ)
    model = generate_series60(p)
    fb_mesh = with_freeboard(model.mesh, freeboard=FREEBOARD, waterline_z=0.0)
    stab = initial_stability(fb_mesh, kg=KG, waterline_z=0.0)
    curve = righting_arm_curve(fb_mesh, kg=KG, heels_deg=HEELS_DEG, waterline_z=0.0)
    
    # Slope at phi=0 should equal GM
    phi1_rad = math.radians(curve.points[1].heel_deg)
    gz1 = curve.points[1].righting_arm
    slope = gz1 / math.sin(phi1_rad)
    assert abs(slope - stab.gm_transverse) / stab.gm_transverse < 0.01, (
        f"Slope {slope:.6f} != GM {stab.gm_transverse:.6f} "
        f"(ratio {slope/stab.gm_transverse:.4f})"
    )


def test_wall_sided_formula_cb070():
    """Wall-sided formula agreement at small angles for Cb=0.70."""
    p = Series60Parameters(length=L, beam=B, draft=T, cb=0.70, nx=NX, nz=NZ)
    model = generate_series60(p)
    fb_mesh = with_freeboard(model.mesh, freeboard=FREEBOARD, waterline_z=0.0)
    stab = initial_stability(fb_mesh, kg=KG, waterline_z=0.0)
    curve = righting_arm_curve(fb_mesh, kg=KG, heels_deg=HEELS_DEG, waterline_z=0.0)
    
    bm = stab.gm_transverse + stab.kb - stab.kg
    # Check at 10 degrees (small angle, wall-sided should be accurate)
    pt_10 = [pt for pt in curve.points if pt.heel_deg == 10][0]
    gz_ws = wall_sided_gz(math.radians(10), stab.gm_transverse, bm)
    assert abs(pt_10.righting_arm - gz_ws) / gz_ws < 0.01, (
        f"GZ {pt_10.righting_arm:.6f} != wall-sided {gz_ws:.6f} "
        f"(error {abs(pt_10.righting_arm - gz_ws) / gz_ws * 100:.2f}%)"
    )
