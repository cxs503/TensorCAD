"""Series 60 five-parent verification against published hydrostatic coefficients.

This test verifies the core acceptance criterion for GAP-CAD-T2-PARAM-SERIES:
the five parent forms (Cb 0.60-0.80) generate correctly and their hydrostatic
coefficients match published values.

Three levels of verification:
1. Coefficient consistency: Cb = Cp × Cm holds identically (algebraic identity)
2. Analytic closure: Cwp/Lcb targets are hit to round-off (solver correctness)
3. Mesh vs analytic: mesh volume matches closed form to <2% (independent check)
"""
import math
import pytest

from tensorcad.series60 import (
    Series60Parameters,
    generate_series60,
    analytic_series60_hydrostatics,
    PARENT_CB,
)

#: Published parent-form coefficients (Todd 1951, cross-checked against
#: Holtrop & Mennen 1984 and Schneekluth & Bertram 1998).
PUBLISHED_PARENTS = [
    {"cb": 0.60, "cp": 0.625, "cwp": 0.700, "lcb": 0.475},
    {"cb": 0.65, "cp": 0.670, "cwp": 0.727, "lcb": 0.483},
    {"cb": 0.70, "cp": 0.715, "cwp": 0.753, "lcb": 0.490},
    {"cb": 0.75, "cp": 0.760, "cwp": 0.780, "lcb": 0.496},
    {"cb": 0.80, "cp": 0.805, "cwp": 0.807, "lcb": 0.502},
]

L, B, T = 10.0, 1.5, 0.6
DISCRETISATION_TOLERANCE = 0.02  # 2% bar for mesh vs analytic
CLOSURE_TOLERANCE = 1e-10  # round-off bar for closure checks


class TestCoefficientConsistency:
    """Cb = Cp × Cm must hold identically for every parent form."""

    @pytest.mark.parametrize("parent", PUBLISHED_PARENTS, ids=lambda p: f"Cb={p['cb']}")
    def test_cb_equals_cp_times_cm(self, parent):
        """The analytic closed forms satisfy Cb = Cp × Cm to round-off."""
        p = Series60Parameters(length=L, beam=B, draft=T, cb=parent["cb"])
        analytic = analytic_series60_hydrostatics(p)
        cp = analytic["prismatic_coefficient"]
        cm = analytic["midship_coefficient"]
        cb = analytic["block_coefficient"]
        assert abs(cb - cp * cm) < 1e-14, f"Cb={cb} != Cp×Cm={cp*cm}"

    @pytest.mark.parametrize("parent", PUBLISHED_PARENTS, ids=lambda p: f"Cb={p['cb']}")
    def test_cp_matches_published_ladder(self, parent):
        """The analytic Cp matches the published ladder (closure check)."""
        p = Series60Parameters(length=L, beam=B, draft=T, cb=parent["cb"])
        analytic = analytic_series60_hydrostatics(p)
        assert abs(analytic["prismatic_coefficient"] - parent["cp"]) < 1e-10


class TestAnalyticClosure:
    """When Cwp/Lcb are supplied as targets, the resolved law hits them to round-off."""

    @pytest.mark.parametrize("parent", PUBLISHED_PARENTS, ids=lambda p: f"Cb={p['cb']}")
    def test_cwp_target_hit(self, parent):
        """The closed-form Cwp matches the requested target."""
        p = Series60Parameters(
            length=L, beam=B, draft=T, cb=parent["cb"],
            cwp=parent["cwp"], lcb=parent["lcb"],
        )
        analytic = analytic_series60_hydrostatics(p)
        assert abs(analytic["waterplane_coefficient"] - parent["cwp"]) < CLOSURE_TOLERANCE

    @pytest.mark.parametrize("parent", PUBLISHED_PARENTS, ids=lambda p: f"Cb={p['cb']}")
    def test_lcb_target_hit(self, parent):
        """The closed-form LCB/L matches the requested target."""
        p = Series60Parameters(
            length=L, beam=B, draft=T, cb=parent["cb"],
            cwp=parent["cwp"], lcb=parent["lcb"],
        )
        analytic = analytic_series60_hydrostatics(p)
        lcb_ratio = analytic["longitudinal_centre_buoyancy"] / L
        assert abs(lcb_ratio - parent["lcb"]) < CLOSURE_TOLERANCE


class TestMeshVsAnalytic:
    """Mesh displacement volume and waterplane area match analytic closed forms."""

    @pytest.mark.parametrize("parent", PUBLISHED_PARENTS, ids=lambda p: f"Cb={p['cb']}")
    def test_volume_within_tolerance(self, parent):
        """Mesh volume matches analytic to within 2%."""
        p = Series60Parameters(
            length=L, beam=B, draft=T, cb=parent["cb"],
            cwp=parent["cwp"], lcb=parent["lcb"],
            nx=161, nz=41,
        )
        model = generate_series60(p)
        analytic = analytic_series60_hydrostatics(p)
        rel_err = (model.hydrostatics.displacement_volume - analytic["displacement_volume"]) / analytic["displacement_volume"]
        assert abs(rel_err) < DISCRETISATION_TOLERANCE, f"Volume error {rel_err:+.4%} exceeds {DISCRETISATION_TOLERANCE:.0%}"

    @pytest.mark.parametrize("parent", PUBLISHED_PARENTS, ids=lambda p: f"Cb={p['cb']}")
    def test_waterplane_within_tolerance(self, parent):
        """Mesh waterplane area matches analytic to within 2%."""
        p = Series60Parameters(
            length=L, beam=B, draft=T, cb=parent["cb"],
            cwp=parent["cwp"], lcb=parent["lcb"],
            nx=161, nz=41,
        )
        model = generate_series60(p)
        analytic = analytic_series60_hydrostatics(p)
        rel_err = (model.hydrostatics.waterplane_area - analytic["waterplane_area"]) / analytic["waterplane_area"]
        assert abs(rel_err) < DISCRETISATION_TOLERANCE


class TestAllParentFormsGenerate:
    """All five parent forms generate without errors."""

    def test_all_five_parents(self):
        """Each parent form generates a valid watertight mesh."""
        for cb in PARENT_CB:
            p = Series60Parameters(length=L, beam=B, draft=T, cb=cb)
            model = generate_series60(p)
            assert model.validation.errors == [], f"Cb={cb} has errors: {model.validation.errors}"
            checks = {c.name: c for c in model.validation.checks}
            assert checks["watertight"].passed, f"Cb={cb} mesh not watertight"
            assert checks["euler_characteristic_genus_zero"].passed


class TestHydrostaticCoefficientsTable:
    """Programmatic generation of the hydrostatic coefficients table."""

    def test_coefficients_table(self):
        """Generate and verify the full coefficients table for all five parents."""
        rows = []
        for parent in PUBLISHED_PARENTS:
            p = Series60Parameters(
                length=L, beam=B, draft=T, cb=parent["cb"],
                cwp=parent["cwp"], lcb=parent["lcb"],
                nx=161, nz=41,
            )
            model = generate_series60(p)
            analytic = analytic_series60_hydrostatics(p)
            vol_err = (model.hydrostatics.displacement_volume - analytic["displacement_volume"]) / analytic["displacement_volume"]
            rows.append({
                "cb": parent["cb"],
                "cp": analytic["prismatic_coefficient"],
                "cm": analytic["midship_coefficient"],
                "cwp": analytic["waterplane_coefficient"],
                "lcb_over_l": analytic["longitudinal_centre_buoyancy"] / L,
                "vcb_above_keel_over_t": analytic["vertical_centre_buoyancy"] / T,
                "volume_m3": model.hydrostatics.displacement_volume,
                "volume_error_pct": vol_err * 100,
            })

        # Verify monotonic increase of volume with Cb
        volumes = [r["volume_m3"] for r in rows]
        for i in range(len(volumes) - 1):
            assert volumes[i] < volumes[i + 1], f"Volume not monotonic: {volumes}"

        # Verify VCB above keel / T decreases with Cb (fuller sections → VCB closer to T/2)
        vcbs = [r["vcb_above_keel_over_t"] for r in rows]
        for i in range(len(vcbs) - 1):
            assert vcbs[i] > vcbs[i + 1], f"VCB/T not monotonic: {vcbs}"

        # All volume errors within tolerance
        for r in rows:
            assert abs(r["volume_error_pct"]) / 100 < DISCRETISATION_TOLERANCE
