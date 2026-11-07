"""Fidelity of the Series 60 family to the published parent forms.

The family was checked against its own closed forms (``tests/test_series60.py``), which
is necessary and not sufficient: a parametrisation can be algebraically exact and still
describe a hull no Series 60 parent resembles.  These tests close that gap.

The finding they pin down: with a single longitudinal exponent the deck width and the
section area are both proportional to ``1 - |2*xi-1|^n``, so ``Cwp = Cp`` and
``LCB = L/2`` *identically*.  Published parents carry ``Cwp`` well above ``Cp`` and their
LCB aft of midship, so the legacy hull is 10.7 % low on ``Cwp`` at ``Cb = 0.60``.
``Series60Parameters.cwp``/``.lcb`` make both quantities inputs; the legacy behaviour is
unchanged when they are omitted.
"""
from __future__ import annotations

import math

import pytest

from tensorcad.series60 import (
    Series60Parameters,
    _resolve_halves,
    analytic_series60_hydrostatics,
    gauss_legendre_unit,
    generate_series60,
    series60_section_law,
)

# Tabulated Series 60 parent-form particulars (Todd 1951, as reproduced in the standard
# resistance references).  ``Cb = Cp * Cm`` holds in the table to the digit shown, which
# is what makes the transcription auditable.
PUBLISHED: dict[float, dict[str, float]] = {
    0.60: {"cm": 0.960, "cp": 0.625, "cwp": 0.700},
    0.65: {"cm": 0.970, "cp": 0.670, "cwp": 0.727},
    0.70: {"cm": 0.980, "cp": 0.715, "cwp": 0.753},
    0.75: {"cm": 0.985, "cp": 0.760, "cwp": 0.780},
    0.80: {"cm": 0.990, "cp": 0.805, "cwp": 0.807},
}
PARENTS = sorted(PUBLISHED)
L, B, T = 10.0, 1.5, 0.6


def schneekluth_cwp(cb: float) -> float:
    """``(1 + 2 Cb) / 3`` -- Schneekluth & Bertram's waterplane regression."""
    return (1.0 + 2.0 * cb) / 3.0


def _quad(fn, order: int = 96) -> float:
    nodes, weights = gauss_legendre_unit(order)
    return math.fsum(w * fn(x) for x, w in zip(nodes, weights))


class TestLegacyStructuralDefect:
    """What a single longitudinal exponent cannot express, proved in-repo."""

    @pytest.mark.parametrize("cb", PARENTS)
    def test_cwp_equals_cp_identically(self, cb: float) -> None:
        ref = analytic_series60_hydrostatics(Series60Parameters(length=L, beam=B, draft=T, cb=cb))
        assert ref["waterplane_coefficient"] == pytest.approx(ref["prismatic_coefficient"], rel=1e-15)

    @pytest.mark.parametrize("cb", PARENTS)
    def test_lcb_is_pinned_to_midship(self, cb: float) -> None:
        law = series60_section_law(Series60Parameters(length=L, beam=B, draft=T, cb=cb))
        assert law.lcb_ratio == pytest.approx(0.5, abs=1e-15)
        assert law.lcf_ratio == pytest.approx(0.5, abs=1e-15)

    def test_cwp_is_far_below_the_published_parents(self) -> None:
        for cb in (0.60, 0.65, 0.70):
            ref = analytic_series60_hydrostatics(Series60Parameters(length=L, beam=B, draft=T, cb=cb))
            error = ref["waterplane_coefficient"] / PUBLISHED[cb]["cwp"] - 1.0
            assert error < -0.03, f"Cb={cb}: unexpected Cwp fidelity ({error:+.2%})"
        worst = analytic_series60_hydrostatics(
            Series60Parameters(length=L, beam=B, draft=T, cb=0.60))["waterplane_coefficient"] / 0.700 - 1.0
        assert worst == pytest.approx(-0.107, abs=0.002)

    def test_cwp_is_far_below_the_regression_at_every_parent(self) -> None:
        for cb in PARENTS:
            ref = analytic_series60_hydrostatics(Series60Parameters(length=L, beam=B, draft=T, cb=cb))
            error = ref["waterplane_coefficient"] / schneekluth_cwp(cb) - 1.0
            assert error < -0.05, f"Cb={cb}: {error:+.2%} against the regression"

    @pytest.mark.parametrize("cb", PARENTS)
    def test_published_table_is_self_consistent(self, cb: float) -> None:
        row = PUBLISHED[cb]
        assert row["cp"] * row["cm"] == pytest.approx(cb, abs=0.004)


class TestClosedFormsAgainstIndependentIntegration:
    """The closed forms are checked by a rule that shares no code with them."""

    @pytest.mark.parametrize("pair", [(1.6667, 1.6667, 0.625), (2.9705, 1.8724, 0.625),
                                      (4.06, 4.31, 0.805), (2.0, 1.4, 0.660)])
    def test_waterplane_and_lcb_match_quadrature(self, pair) -> None:
        n_aft, n_fore, cp = pair

        def half(fn, n: float) -> float:
            return _quad(lambda t: fn((1.0 - t) / 2.0, n))

        def fore(fn, n: float) -> float:
            return _quad(lambda t: fn((1.0 + t) / 2.0, n))

        def taper(xi, n):
            return 1.0 - abs(2.0 * xi - 1.0) ** n

        cwp_hatch = 0.5 * (half(lambda xi, n: taper(xi, n), n_aft)
                           + fore(lambda xi, n: taper(xi, n), n_fore))
        square = 0.5 * (half(lambda xi, n: taper(xi, n) ** 2, n_aft)
                        + fore(lambda xi, n: taper(xi, n) ** 2, n_fore))
        first = 0.5 * (half(lambda xi, n: xi * taper(xi, n), n_aft)
                       + fore(lambda xi, n: xi * taper(xi, n), n_fore))
        first_square = 0.5 * (half(lambda xi, n: xi * taper(xi, n) ** 2, n_aft)
                              + fore(lambda xi, n: xi * taper(xi, n) ** 2, n_fore))
        gain = (cp - cwp_hatch) / (cwp_hatch - square)
        lcb_hatch = ((1.0 + gain) * first - gain * first_square) / cp

        cwp, lcb, _square, _gain = _resolve_halves(n_aft, n_fore, cp)
        # 96-node Gauss-Legendre on integrands carrying t**1.4 and t**4.3: agreement
        # is limited by the rule, not by the closed forms, hence the loose-looking tolerance.
        assert cwp == pytest.approx(cwp_hatch, rel=1e-9)
        assert lcb == pytest.approx(lcb_hatch, rel=1e-8)

    def test_lcf_is_the_waterplane_centroid(self) -> None:
        for cb, cwp, lcb in ((0.60, 0.700, 0.475), (0.75, 0.780, 0.496)):
            law = series60_section_law(Series60Parameters(length=L, beam=B, draft=T, cb=cb,
                                                          cwp=cwp, lcb=lcb))
            def taper(xi, n):
                return 1.0 - abs(2.0 * xi - 1.0) ** n
            moment = _quad(lambda t: 0.5 * ((1 - t) / 2 * taper((1 - t) / 2, law.n_aft)
                                            + (1 + t) / 2 * taper((1 + t) / 2, law.n_fore)))
            area = _quad(lambda t: 0.5 * (taper((1 - t) / 2, law.n_aft)
                                          + taper((1 + t) / 2, law.n_fore)))
            assert law.lcf_ratio == pytest.approx(moment / area, rel=1e-10)

    @pytest.mark.parametrize("kwargs", [{}, {"cwp": 0.700}, {"cwp": 0.700, "lcb": 0.475}])
    def test_vcb_matches_independent_vertical_quadrature(self, kwargs) -> None:
        law = series60_section_law(Series60Parameters(length=L, beam=B, draft=T, cb=0.60, **kwargs))

        def station_value(xi):
            q = law.fullness(xi)
            m = q / (1.0 - q)
            # int_0^1 zeta**(1/m) dzeta = m/(m+1); int_0^1 zeta*zeta**(1/m) dzeta = m/(2m+1)
            return m

        def vcb_density(t):
            for side, xi in ((law.n_aft, (1 - t) / 2), (law.n_fore, (1 + t) / 2)):
                pass
            return None

        # Volume-weighted mean of VCB/T, integrated over t for both halves.
        def one_half(t, n):
            xi = (1 - t) / 2
            f = 1.0 - t ** n
            q = law.cm_mid * (1.0 + law.gain * (1.0 - f))
            m = q / (1.0 - q)
            return f * q * ((m + 1.0) / (2.0 * m + 1.0))

        def other_half(t, n):
            xi = (1 + t) / 2
            f = 1.0 - t ** n
            q = law.cm_mid * (1.0 + law.gain * (1.0 - f))
            m = q / (1.0 - q)
            return f * q * ((m + 1.0) / (2.0 * m + 1.0))

        hatch = 0.5 * (_quad(lambda t: one_half(t, law.n_aft))
                       + _quad(lambda t: other_half(t, law.n_fore)))
        assert law.vcb_ratio == pytest.approx(hatch / law.cb, rel=1e-10)


class TestIndependentParticulars:
    """Cwp and LCB become inputs, and the loft follows them."""

    @pytest.mark.parametrize("cb", PARENTS)
    def test_published_particulars_are_reachable(self, cb: float) -> None:
        target = PUBLISHED[cb]["cwp"]
        lcb_target = 0.5 - 0.025 if cb == 0.60 else 0.5
        p = Series60Parameters(length=L, beam=B, draft=T, cb=cb, cwp=target, lcb=lcb_target)
        law = series60_section_law(p)
        ref = analytic_series60_hydrostatics(p)
        assert ref["waterplane_coefficient"] == pytest.approx(target, rel=1e-13)
        assert law.lcb_ratio == pytest.approx(lcb_target, abs=1e-13)
        assert ref["midship_coefficient"] == pytest.approx(PUBLISHED[cb]["cm"], abs=0.01)
        assert ref["block_coefficient"] == pytest.approx(cb, rel=1e-13)

    def test_gain_zero_reproduces_the_legacy_hull(self) -> None:
        legacy = Series60Parameters(length=L, beam=B, draft=T, cb=0.60, nx=41, nz=13)
        implied = analytic_series60_hydrostatics(legacy)["waterplane_coefficient"]
        explicit = Series60Parameters(length=L, beam=B, draft=T, cb=0.60, nx=41, nz=13,
                                      cwp=implied, lcb=0.5)
        assert series60_section_law(explicit).gain == pytest.approx(0.0, abs=1e-12)
        a = generate_series60(legacy)
        b = generate_series60(explicit)
        for va, vb in zip(a.mesh.vertices, b.mesh.vertices):
            for ca, cb_ in zip(va, vb):
                assert cb_ == pytest.approx(ca, rel=1e-12, abs=1e-15)

    def test_mesh_reproduces_the_requested_particulars(self) -> None:
        p = Series60Parameters(length=L, beam=B, draft=T, cb=0.60, nx=161, nz=41,
                               cwp=0.700, lcb=0.475)
        model = generate_series60(p)
        ref = analytic_series60_hydrostatics(p)
        hydro = model.hydrostatics
        assert hydro.waterplane_area / ref["waterplane_area"] - 1.0 == pytest.approx(0.0, abs=1e-3)
        assert hydro.displacement_volume / ref["displacement_volume"] - 1.0 == pytest.approx(0.0, abs=0.02)
        assert (hydro.longitudinal_centre_buoyancy / L - 0.475) == pytest.approx(0.0, abs=2e-4)
        assert model.validation.errors == []

    def test_mesh_lcb_error_shrinks_under_refinement(self) -> None:
        errors = []
        for nx, nz in ((81, 21), (161, 41)):
            p = Series60Parameters(length=L, beam=B, draft=T, cb=0.60, nx=nx, nz=nz,
                                   cwp=0.700, lcb=0.475)
            model = generate_series60(p)
            errors.append(abs(model.hydrostatics.longitudinal_centre_buoyancy / L - 0.475))
        assert errors[1] < 0.6 * errors[0], errors

    def test_controlled_hull_stays_a_closed_manifold(self) -> None:
        p = Series60Parameters(length=L, beam=B, draft=T, cb=0.65, nx=81, nz=21,
                               cwp=0.727, lcb=0.483)
        model = generate_series60(p)
        assert model.validation.errors == []
        for check in model.validation.checks:
            if check.name in ("watertight", "euler_characteristic_genus_zero",
                              "outward_orientation"):
                assert check.passed, check.message


class TestRejections:
    def test_lcb_without_cwp_is_rejected(self) -> None:
        with pytest.raises(ValueError, match="lcb needs cwp"):
            Series60Parameters(cb=0.60, lcb=0.475).validate()

    def test_waterplane_too_fine_for_the_block_coefficient_is_rejected(self) -> None:
        with pytest.raises(ValueError, match="section coefficient"):
            Series60Parameters(cb=0.70, cwp=0.700).validate()

    def test_out_of_range_cwp_is_rejected(self) -> None:
        with pytest.raises(ValueError):
            Series60Parameters(cb=0.60, cwp=0.30).validate()

    def test_gauss_legendre_integrates_polynomials_exactly(self) -> None:
        nodes, weights = gauss_legendre_unit(8)
        for power in range(2 * len(nodes)):
            exact = 1.0 / (power + 1)
            got = math.fsum(w * x ** power for x, w in zip(nodes, weights))
            assert got == pytest.approx(exact, rel=1e-13, abs=1e-15)