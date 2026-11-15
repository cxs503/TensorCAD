"""Series 60: the declared Cb window, and a reference for the wetted surface area.

Two separate debts are closed here.

**The declared domain.**  ``Series60Parameters`` and ``Series60Parameters.problems``
both advertised ``Cb`` in ``[0.55, 0.85]``, but the published Todd ladder ends at
``Cb = 0.80`` with ``Cp = 0.805``, and ``Cm = Cb / Cp`` needs a prismatic coefficient
above the block coefficient.  Every request in ``(0.805, 0.85]`` therefore raised
``Cp (0.805000) must exceed Cb (0.850000)`` from *inside* the interpolation -- an
advertised, unreachable corner of the parameter space.  The window is now the ladder
(``CB_MIN``/``CB_MAX``) and this module sweeps the whole of it, boundary included,
because that is what the advertising claim means.

**The wetted surface.**  ``wetted_surface_area`` is the one quantity the family had no
reference for at all (``NO_CLOSED_FORM``), and it is the quantity the resistance chain
consumes.  :func:`tensorcad.series60.series60_surface_area_reference` integrates the
parametric metric two independent ways -- 2-D Gauss-Legendre on the analytic
Jacobian, and parametric quads on a refinement ladder, Richardson-extrapolated -- and
this module checks the two against each other, against the lofted mesh at the *same*
grid, and against the mesh under refinement.
"""
import math

import pytest

from tensorcad.series60 import (
    CB_MAX,
    CB_MIN,
    CP_MAX,
    PARENT_CB,
    Series60Parameters,
    _deck_taper_derivative,
    _exponent_derivative,
    build_series60_mesh,
    hull_surface_area_quad_sum,
    hull_surface_area_quadrature,
    series60_section_law,
    series60_surface_area_audit,
    series60_surface_area_reference,
)
from tensorcad.hydrostatics import compute_hydrostatics

CASE = dict(length=10.0, beam=1.5, draft=0.6)


def _model(cb=0.70, nx=161, nz=41, **extra):
    return Series60Parameters(cb=cb, nx=nx, nz=nz, **CASE, **extra)


# --------------------------------------------------------------------------------------
# The declared Cb window is the published ladder, and every point of it generates.
# --------------------------------------------------------------------------------------

def test_the_published_ladder_is_inside_the_declared_window():
    assert CB_MIN <= min(PARENT_CB) and max(PARENT_CB) <= CB_MAX
    assert CB_MAX == max(PARENT_CB)


def test_every_declared_block_coefficient_generates_and_hits_its_target():
    """The property the old window failed: advertised must imply reachable.

    Sweeps ``[CB_MIN, CB_MAX]`` inclusive at 0.005 steps, i.e. both boundaries and
    every interpolated parent.  A coarse grid is deliberate -- this is a domain test,
    not a discretisation test, so the resolved *law* is the thing compared.
    """
    cb = CB_MIN
    checked = []
    while cb <= CB_MAX + 1e-12:
        law = series60_section_law(_model(cb=round(cb, 6), nx=21, nz=7))
        assert law.cb == pytest.approx(round(cb, 6), rel=0, abs=1e-9), f"Cb={cb}"
        assert law.cwp > law.cb, f"Cwp must exceed Cb at Cb={cb}"
        assert law.cm_mid <= 1.0 and law.cm_mid > law.cb, f"Cm=Cb/Cp ill-formed at Cb={cb}"
        checked.append(cb)
        cb += 0.005
    assert len(checked) > 45
    assert checked[0] == pytest.approx(CB_MIN) and checked[-1] == pytest.approx(CB_MAX, abs=1e-9)


def test_the_old_dead_zone_now_says_why_instead_of_leaking_the_cp_check():
    """``Cb`` above the ladder is refused by name, not by an internal assertion.

    Regression for the defect: at ``Cb = 0.82`` the message used to be
    ``Cp (0.805000) must exceed Cb (0.820000)``, which names neither the parameter
    nor its admissible range.
    """
    for cb in (0.801, 0.82, 0.85, 0.9, 0.4):
        with pytest.raises(ValueError) as caught:
            series60_section_law(_model(cb=cb))
        message = str(caught.value)
        assert "Cb must lie in" in message and str(CB_MAX) in message, message
        assert "Todd" in message, message
    # and the same through the parameter object's own validator
    problems, _warnings = Series60Parameters(cb=0.85, **CASE).problems()
    assert problems and "cb must lie in" in problems[0].lower()


def test_the_boundary_of_the_window_is_the_fullest_published_parent():
    """``Cb = CB_MAX`` needs ``Cp > Cb``; ``CP_MAX`` is what makes that reachable."""
    assert CP_MAX > CB_MAX
    law = series60_section_law(_model(cb=CB_MAX, nx=21, nz=7))
    assert law.cm_mid == pytest.approx(CB_MAX / CP_MAX, rel=1e-12)


# --------------------------------------------------------------------------------------
# The surface reference: its derivatives, its two paths, and the loft.
# --------------------------------------------------------------------------------------

@pytest.mark.parametrize("cb,transom", [(0.65, 0.0), (0.75, 0.15), (0.80, 0.25)])
def test_analytic_derivatives_agree_with_finite_differences(cb, transom):
    """The metric needs ``f'`` and ``m'``; both are hand-derived, so both are tested."""
    law = series60_section_law(_model(cb=cb, transom_ratio=transom))
    h = 1e-6
    for xi in (0.05, 0.2, 0.35, 0.6, 0.75, 0.95):
        numeric_f = (law.deck_taper(xi + h) - law.deck_taper(xi - h)) / (2.0 * h)
        assert _deck_taper_derivative(law, xi) == pytest.approx(numeric_f, rel=1e-7, abs=1e-9)
        numeric_m = (law.exponent(xi + h) - law.exponent(xi - h)) / (2.0 * h)
        assert _exponent_derivative(law, xi) == pytest.approx(numeric_m, rel=1e-7, abs=1e-6)


@pytest.mark.parametrize("cb,transom", [(0.60, 0.0), (0.70, 0.15), (0.80, 0.0)])
def test_the_two_reference_paths_agree(cb, transom):
    """Quadrature (analytic Jacobian) vs extrapolated quad sum (corner areas).

    They share the station law and nothing else, so agreement is evidence about the
    reference; disagreement would mean the reference is wrong, not the mesh.
    """
    reference = series60_surface_area_reference(_model(cb=cb, transom_ratio=transom))
    assert reference["path_agreement"] < 5e-4, reference
    # the ks near 170 (the fullest parent) make the rule converge slowly: 1e-4 is
    # the measured self-error bar, and the second path is what actually pins the value
    assert reference["hull_area_quadrature_error_m2"] < 1e-4 * reference["hull_area_extrapolated_m2"]
    assert reference["area_quad_sum_agreement"] < 1e-2


def test_the_reference_integrates_the_same_surface_the_loft_polygonises():
    """At the *same* grid, parametric quads and the mesh triangles must coincide.

    This is the test that ties the reference to production: an error in the metric
    (a wrong ``f'`` sign, the wrong vertical exponent, a missing transom face) shows
    up here as a difference at equal resolution, before any convergence question.
    """
    p = _model(cb=0.70, nx=81, nz=31, transom_ratio=0.12)
    mesh, _offsets = build_series60_mesh(p)
    hydro = compute_hydrostatics(mesh, 0.0, water_density=p.water_density, closed=True)
    law = series60_section_law(p)
    # The mesh grid samples stations 0..nx-1 and levels 0..nz-1 of the *same* surface.
    k = 1.0
    from tensorcad.series60 import vertical_exponent_range
    k = vertical_exponent_range(p, law)[1] * 1.05 + 1e-9
    # levels are zeta = level/(nz-1), so u = zeta**(1/k) reproduces them exactly
    nz, nx = p.nz, p.nx

    def u_of(level):
        return (level / (nz - 1)) ** (1.0 / k)

    from tensorcad.series60 import _surface_point, _triangle_area
    total = 0.0
    rows = [[_surface_point(p, law, i / (nx - 1), u_of(level), k) for level in range(nz)]
            for i in range(nx)]
    for i in range(nx - 1):
        for level in range(nz - 1):
            a, b = rows[i][level], rows[i + 1][level]
            c, d = rows[i + 1][level + 1], rows[i][level + 1]
            total += 2.0 * (_triangle_area(a, b, c) + _triangle_area(a, c, d))
    transom = hydro.wetted_surface_area
    # the mesh includes the stern face; the parametric sum above does not
    reference_face = series60_surface_area_reference(p)["transom_area_m2"]
    assert total + reference_face == pytest.approx(transom, rel=2e-2)


@pytest.mark.parametrize("cb", [0.60, 0.70, 0.80])
def test_mesh_wetted_area_converges_to_the_reference(cb):
    """The measured ladder closes on the reference from below, and it must.

    A polygonised hull is inscribed in the surface, so its area is short of the true
    area; the bar is the same 2 % the module uses for volume, applied at the grid the
    CLI defaults to (161 x 41).
    """
    reference = series60_surface_area_reference(_model(cb=cb))["total_reference_m2"]
    errors = {}
    for nx, nz in ((41, 13), (81, 21), (161, 41), (321, 81)):
        audit = series60_surface_area_audit(_model(cb=cb, nx=nx, nz=nz))
        errors[(nx, nz)] = audit["relative_error"]
        assert audit["measured_m2"] < reference * (1.0 + 1e-12), (nx, nz)
    assert errors[(161, 41)] > -0.02, errors
    assert abs(errors[(321, 81)]) < abs(errors[(161, 41)]) < abs(errors[(41, 13)]), errors


def test_a_transom_stays_inside_the_surface_bar():
    """The ends are where a wetted-surface bookkeeping error would hide."""
    plain = series60_surface_area_audit(_model(cb=0.70, nx=161, nz=41))
    transom = series60_surface_area_audit(_model(cb=0.70, nx=161, nz=41, transom_ratio=0.20))
    for audit in (plain, transom):
        assert abs(audit["relative_error"]) < 0.02, audit
    # the stern face is genuinely part of the reference, not a constant
    assert transom["reference"]["transom_area_m2"] > 0.0
    assert transom["reference_m2"] > plain["reference_m2"]
    assert transom["reference"]["transom_area_m2"] == pytest.approx(
        transom["reference_m2"] - transom["reference"]["hull_area_extrapolated_m2"], rel=1e-12
    )