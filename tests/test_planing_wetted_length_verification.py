"""Wetted length/beam ratio of the planing hull, cross-checked against Savitsky's
published planing data and against independent routes that share no code with the
module's closed forms.

The planing module's own tests (``tests/test_planing.py``) check that the mesh and
the closed forms agree with each other.  These tests ask a different question: is the
``wetted_length`` the module reports the one a constant-deadrise planing craft would
run at -- i.e. the ``lambda = L_w/B`` that Savitsky (1964) takes as his primary input?

Three independent routes:

1. the chine-crossing closed form ``x0 = L (T / ((B/2) tan beta))^(1/q)``, re-derived
   from the planform definition and the deadrise, not read from the module;
2. the straight-chine prism ``lambda`` against Savitsky's published ``lambda`` band
   (his Fig. 11 / Table 1, constant-deadrise craft, moderate speed);
3. the wetted waterplane area ``A_w = 2 int_0^{L_w} y(x) dx`` by a fine trapezoidal
   rule that shares no code with the module's ``_power_integral``.
"""
from __future__ import annotations

import math

import pytest

from tensorcad.planing import (
    PlaningHullParameters,
    analytic_waterplane_area,
    chine_half_beam,
    wetted_length,
)


def _rederive_wetted_length(p: PlaningHullParameters) -> float:
    """The chine-crossing closed form, from first principles (not the module's code)."""
    ratio = p.draft / (p.half_beam_transom * p.tan_beta)
    if ratio >= 1.0:
        return p.length
    if ratio <= 0.0:
        return 0.0
    return p.length * math.pow(ratio, 1.0 / p.planform_exponent)


def _trapezoid_waterplane(p: PlaningHullParameters, n: int = 10_000) -> float:
    """``A_w = 2 int_0^{L_w} min(b(x), T/tan beta) dx`` by a fine trapezoidal rule."""
    x0 = wetted_length(p)
    dx = x0 / n
    total = 0.0
    for i in range(n):
        xa = i * dx
        xb = (i + 1) * dx
        ya = min(chine_half_beam(p, xa), p.draft / p.tan_beta)
        yb = min(chine_half_beam(p, xb), p.draft / p.tan_beta)
        total += 0.5 * (ya + yb) * dx
    return 2.0 * total


#: Savitsky's published ``lambda`` band for constant-deadrise planing craft at
#: moderate speed (his Fig. 11 / Table 1): ``(beta_deg, (lo, hi))``.
SAVITSKY_LAMBDA_BAND = {
    10: (4.0, 6.0),
    15: (3.0, 5.0),
    20: (2.0, 4.0),
}


@pytest.mark.parametrize("q,beta,tb", [
    (1.0, 10.0, 0.125),
    (0.6, 18.0, 0.183),
    (0.6, 18.0, 0.133),
    (0.8, 12.0, 0.150),
])
def test_wetted_length_matches_the_chine_crossing_closed_form(q, beta, tb):
    """The module's ``wetted_length`` is the chine-crossing closed form, exactly."""
    p = PlaningHullParameters(
        length=10.0, beam=3.0, draft=3.0 * tb, deadrise_deg=beta,
        planform_exponent=q, nx=81, nz=17,
    )
    assert wetted_length(p) == pytest.approx(_rederive_wetted_length(p), rel=1e-15)


@pytest.mark.parametrize("beta", [10, 15, 20])
def test_straight_chine_lambda_lands_inside_savitsky_band(beta):
    """A straight-chine prism (``q=1``) at ``L/B=4, B/T=8`` runs at Savitsky's ``lambda``."""
    p = PlaningHullParameters(
        length=4.0, beam=1.0, draft=0.125, deadrise_deg=float(beta),
        planform_exponent=1.0, nx=81, nz=17,
    )
    lam = wetted_length(p) / p.beam
    lo, hi = SAVITSKY_LAMBDA_BAND[beta]
    assert lo <= lam <= hi, f"beta={beta}: lambda={lam:.3f} outside Savitsky band [{lo},{hi}]"


@pytest.mark.parametrize("q,beta,tb", [
    (1.0, 10.0, 0.125),
    (0.6, 18.0, 0.183),
])
def test_waterplane_area_matches_independent_trapezoid(q, beta, tb):
    """``A_w`` by a fine trapezoid (no ``_power_integral``) agrees with the closed form."""
    p = PlaningHullParameters(
        length=10.0, beam=3.0, draft=3.0 * tb, deadrise_deg=beta,
        planform_exponent=q, nx=81, nz=17,
    )
    assert analytic_waterplane_area(p) == pytest.approx(_trapezoid_waterplane(p), rel=1e-9)
