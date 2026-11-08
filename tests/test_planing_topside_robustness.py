"""Regression tests for the planing wetted-topside integral.

Before this test file the integral was

    area = 2 int_0^x0 (T - b tan(beta)) sqrt(1 + b'(x)**2) dx

evaluated by a fixed 4000-interval Simpson rule in the substituted variable
``x = x0 t**(2/q)``, with the square root formed as ``sqrt(1 + slope**2) * dx/dt``.
Two failure modes of that form are reproduced here so they cannot come back:

* **overflow** -- for a fine bowed planing bottom (``q = 0.05``) the chine slope near
  the stem reaches ``~1e21``; squaring it overflows to ``inf`` and the whole wetted
  area becomes ``inf`` on 4 of the 420 parameter sets in
  ``tools/planing_topside_robustness.py`` (``deadrise >= 30 deg``, ``draft <= 0.2 m``).
  An ``inf`` reference silently poisons every claim validated against it.
* **truncation** -- with ``q >= 2`` the ``t**(2/q - 1)`` substitution singularity sits at
  ``t = 0`` and the uniform grid is 0.79-0.85 % short on the 126 such cases; with
  ``q = 0.05`` the arc element is a ``t**39`` boundary layer against ``t = 1``.

The reference values below are ``mpmath.quad`` at 50 digits in a *physical* variable
(``u = b(x)`` for ``q < 1``, ``x`` for ``q >= 1``), cross-checked against
``scipy.integrate.quad`` to ``3e-10`` relative -- see
``docs/planing_numerical_robustness.md``.  The tests are standard-library only: the
numbers are pinned, mpmath is not a test dependency.
"""
from __future__ import annotations

import json
import math

import pytest

from tensorcad import planing
from tensorcad.planing import (
    PlaningHullParameters,
    _hydrostatic_checks,
    analytic_planing_hydrostatics,
    analytic_waterplane_area,
    analytic_wetted_area,
    chine_half_beam,
    generate_planing_hull,
    section_area,
    wetted_length,
)

# (deadrise [deg], planform exponent q, draft [m]) -> wetted area [m^2] at 50 digits.
REFERENCES = {
    # formerly inf: slope**2 overflow (the defect)
    (30.0, 0.05, 0.05): 2.0086602540378447559,
    (45.0, 0.05, 0.05): 1.4192135623730952854,
    (59.0, 0.05, 0.2): 4.7146024383835264038,
    # formerly 0.79-0.85 % low: uniform grid on a t**(2/q - 1) singularity
    (45.0, 4.0, 2.0): 46.753709203666148631,
    (18.0, 4.0, 2.0): 50.46090479207479678,
    (59.0, 3.0, 0.55): 12.191131807710858435,
    # formerly converged: the fix must not move these (docs/planing_hull.md values)
    (18.0, 0.6, 0.55): 25.670434201210412428,
    (18.0, 0.6, 0.4): 21.636391220055819854,
    (1.0, 0.1, 0.35): 35.442136476707574362,
    (12.0, 1.5, 1.0): 32.444512567097557364,
}

SWEEP_DEGREES = (1.0, 18.0, 45.0, 59.0)
SWEEP_EXPONENTS = (0.05, 0.1, 0.6, 1.0, 2.0, 4.0)
SWEEP_DRAFTS = (0.05, 0.55, 2.0)


def parameters(deadrise: float, exponent: float, draft: float, **kw) -> PlaningHullParameters:
    return PlaningHullParameters(
        length=kw.pop("length", 10.0),
        beam=kw.pop("beam", 3.0),
        draft=draft,
        deadrise_deg=deadrise,
        planform_exponent=exponent,
        **kw,
    )


def legacy_topside_area(p: PlaningHullParameters, x0: float, samples: int = 4000) -> float:
    """The pre-0.24 integrand, verbatim, as the "before" side of the regression."""
    if x0 <= 0.0:
        return 0.0
    if samples % 2:
        samples += 1
    exponent = 2.0 / p.planform_exponent
    c = p.half_beam_transom / p.length**p.planform_exponent

    def integrand(t: float) -> float:
        x = x0 * t**exponent
        if x <= 0.0:
            return 0.0
        b = c * x**p.planform_exponent
        slope = c * p.planform_exponent * x ** (p.planform_exponent - 1.0)
        height = p.draft - b * p.tan_beta
        if height <= 0.0:
            return 0.0
        dx_dt = x0 * exponent * t ** (exponent - 1.0)
        return 2.0 * height * math.sqrt(1.0 + slope * slope) * dx_dt

    h = 1.0 / samples
    total = integrand(0.0) + integrand(1.0)
    for i in range(1, samples):
        total += (4.0 if i % 2 else 2.0) * integrand(i * h)
    return total * h / 3.0


def legacy_wetted_area(p: PlaningHullParameters) -> float:
    tau = p.tan_beta
    x0 = wetted_length(p)
    bottom = analytic_waterplane_area(p) * math.sqrt(1.0 + tau * tau)
    return bottom + legacy_topside_area(p, x0) + section_area(p, p.length)


@pytest.mark.parametrize("case,reference", sorted(REFERENCES.items()))
def test_wetted_area_matches_high_precision_reference(case, reference) -> None:
    p = parameters(*case)
    p.validate()
    measured = analytic_wetted_area(p)
    relative = (measured - reference) / reference
    assert abs(relative) < 1e-12, (case, measured, reference, relative)


@pytest.mark.parametrize("case", [(30.0, 0.05, 0.05), (45.0, 0.05, 0.05), (59.0, 0.05, 0.2)])
def test_legacy_overflow_is_reproduced(case) -> None:
    """The defect, pinned: the old integrand returns a non-finite area on these hulls."""
    p = parameters(*case)
    p.validate()
    assert not math.isfinite(legacy_wetted_area(p))
    assert math.isfinite(analytic_wetted_area(p))


@pytest.mark.parametrize("case", [(45.0, 4.0, 2.0), (18.0, 4.0, 2.0), (59.0, 3.0, 0.55)])
def test_legacy_truncation_is_reproduced(case) -> None:
    """The second defect, pinned: the old fixed grid sat 0.1-0.85 % below the truth."""
    p = parameters(*case)
    p.validate()
    reference = REFERENCES[case]
    legacy_error = (legacy_wetted_area(p) - reference) / reference
    assert legacy_error < -1e-4, legacy_error
    assert abs((analytic_wetted_area(p) - reference) / reference) < 1e-12


@pytest.mark.parametrize("deadrise", SWEEP_DEGREES)
@pytest.mark.parametrize("exponent", SWEEP_EXPONENTS)
@pytest.mark.parametrize("draft", SWEEP_DRAFTS)
def test_closed_forms_stay_finite_over_the_family(deadrise, exponent, draft) -> None:
    p = parameters(deadrise, exponent, draft)
    errors, _ = p.problems()
    assert not errors, errors
    reference = analytic_planing_hydrostatics(p)
    for name, value in reference.items():
        assert math.isfinite(value), (name, value)
    for name in (
        "displacement_volume",
        "waterplane_area",
        "wetted_surface_area",
        "transom_section_area",
    ):
        assert reference[name] > 0.0, (name, reference[name])
    assert reference["wetted_surface_area"] > reference["waterplane_area"]


def test_quadrature_tolerance_is_honoured() -> None:
    """A loose tolerance is cheaper but still inside its own band; the default is tight."""
    p = parameters(45.0, 0.05, 2.0)
    tight = analytic_wetted_area(p)
    loose = analytic_wetted_area(p, 1e-6)
    assert math.isfinite(loose)
    assert abs((loose - tight) / tight) < 1e-6
    with pytest.raises(ValueError):
        analytic_wetted_area(p, 0.0)


def test_samples_keyword_still_accepted() -> None:
    """``samples`` was the fixed Simpson count; it is ignored, not an error."""
    p = parameters(18.0, 0.6, 0.55)
    assert analytic_wetted_area(p, samples=200) == pytest.approx(analytic_wetted_area(p))


def test_extreme_bow_flare_generates_and_reports_finite_metadata() -> None:
    """The formerly-``inf`` hull now lofts, validates and serialises strictly."""
    p = parameters(30.0, 0.05, 0.05, nx=41, nz=17)
    model = generate_planing_hull(p)
    assert not model.validation.errors
    reference = analytic_planing_hydrostatics(p)
    assert math.isfinite(reference["wetted_surface_area"])
    measured = model.hydrostatics.wetted_surface_area
    assert abs((measured - reference["wetted_surface_area"]) / reference["wetted_surface_area"]) < 0.02
    # ``allow_nan=False`` is the strict JSON encoder: an ``inf`` reference would raise.
    json.dumps(model.metadata(), allow_nan=False)


def test_broken_analytic_reference_is_an_error_not_a_warning(monkeypatch) -> None:
    """A non-finite closed form can no longer be reported as a soft sampling warning."""
    p = parameters(18.0, 0.6, 0.55, nx=21, nz=9)
    model = generate_planing_hull(p)
    reference = analytic_planing_hydrostatics(p)

    def broken(parameters_):  # noqa: ANN001 - deliberate stand-in
        return {**reference, "wetted_surface_area": float("inf")}

    monkeypatch.setattr(planing, "analytic_planing_hydrostatics", broken)
    checks = _hydrostatic_checks(p, model.hydrostatics)
    wetted = [c for c in checks if c.name == "planing_wetted_area"]
    assert len(wetted) == 1
    assert not wetted[0].passed
    assert wetted[0].severity == "error"
    with pytest.raises(ValueError, match="failed validation"):
        generate_planing_hull(p)


def test_chine_half_beam_is_the_power_law_the_integrand_assumes() -> None:
    """``b(x0 t**e) = b(x0) t**(e q)`` is what makes the arc element overflow-free."""
    p = parameters(30.0, 0.05, 0.05)
    x0 = wetted_length(p)
    b_x0 = chine_half_beam(p, x0)
    for t in (0.1, 0.37, 0.9):
        expected = b_x0 * t ** (2.0 / p.planform_exponent * p.planform_exponent)
        assert chine_half_beam(p, x0 * t ** (2.0 / p.planform_exponent)) == pytest.approx(
            expected, rel=1e-13
        )