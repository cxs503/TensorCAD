"""Parameter-box conformance gate for the Series 60 family.

The family's other test modules pin individual claims: ``tests/test_series60.py`` the
closed forms, ``test_series60_fidelity.py`` the published particulars,
``test_series60_bulb_transom.py`` the ends, ``test_series60_cli.py`` the entry.  What
was missing is the question *does the reported comparison stay well formed across the
whole declared box* -- ``Cb`` 0.60-0.80 x flat transom x bulb x ``(cwp, lcb)``
targets -- because a comparison table that quietly returns ``nan``, or compares a
quantity in the wrong frame, passes every one of those narrower tests.

That gap was not hypothetical.  ``series60_analytic_relative_errors`` used to report
``vertical_centre_buoyancy`` against the closed form measured **up from the keel**
while the mesh integrates ``z`` from the **design waterline**, so every Series 60 hull
in the family came out at about ``-193 %`` on ``VCB``; and it read
``midship_coefficient`` / ``waterplane_coefficient`` off attributes
:class:`~tensorcad.hydrostatics.Hydrostatics` does not carry, so both came out
``nan``.  Nothing failed, because nothing swept the box.  This module is that sweep,
and it is the guard the mutants ``S60-19`` .. ``S60-21`` in
``tools/hull_gate_power/mutate.py`` are run against.

Predeclared thresholds (stated before the numbers, and they are *screens* on a coarse
grid plus one convergence subset, not a substitute for the refinement study in
``docs/series60_verification.md``):

* ``SCREEN_GRID = (81, 21)`` -- wide sweep, bar ``SCREEN_TOLERANCE = 3 %``.  This grid
  is deliberately coarse: the volume error is about ``-2.3 %`` here, so the bar checks
  that the table is well formed and the frames agree, not that it has converged.
* ``FINE_GRID = (161, 41)`` -- the five parents only, bar ``FINE_TOLERANCE = 1.5 %``,
  which is half the module's own ``DISCRETISATION_TOLERANCE = 2 %``.
* ``VCB_VOLUME_COUPLING = 0.5 %`` -- ``VCB`` is a volume-normalised first moment, so
  its discretisation error cannot exceed the volume error by more than half a
  percentage point.  A frame mismatch breaks this by two orders of magnitude.
* ``CLOSURE_TOLERANCE = 0.5 %`` -- ``(cwp, lcb)`` are generator *targets*, solved for,
  so the mesh must reproduce them at screen resolution.
"""
from __future__ import annotations

import math

import pytest

from tensorcad.series60 import (
    PARENT_CB,
    Series60Parameters,
    analytic_series60_hydrostatics,
    generate_series60,
    series60_analytic_relative_errors,
)

L, B, T = 10.0, 1.5, 0.6

SCREEN_GRID = (81, 21)
FINE_GRID = (161, 41)
SCREEN_TOLERANCE = 0.03
FINE_TOLERANCE = 0.015
VCB_VOLUME_COUPLING = 0.005
CLOSURE_TOLERANCE = 0.005

#: Flat-transom ratios spanning the declared range ``0 <= transom_ratio < 0.45``.
TRANSOMS = (0.0, 0.2, 0.4)

#: Published Series 60 waterplane coefficients of the five parents (Todd 1951, the
#: same transcription ``tests/test_series60_fidelity.py`` audits digit by digit).
PUBLISHED_CWP = {0.60: 0.700, 0.65: 0.727, 0.70: 0.753, 0.75: 0.780, 0.80: 0.807}

#: Design longitudinal centre of buoyancy of the family, aft of midship as published.
PUBLISHED_LCB = 0.490

#: Quantities whose measured value must always be compared: the ones a closed-surface
#: integration reports directly.
ALWAYS_COMPARED = (
    "displacement_volume",
    "waterplane_area",
    "maximum_section_area",
    "vertical_centre_buoyancy",
    "block_coefficient",
    "prismatic_coefficient",
    "midship_coefficient",
    "waterplane_coefficient",
)


def _model(cb: float, *, transom: float = 0.0, grid: tuple[int, int] = SCREEN_GRID,
           cwp: float | None = None, lcb: float | None = None,
           bulb: tuple[float, float] | None = None) -> object:
    bulb_length, bulb_area = bulb if bulb is not None else (0.0, 0.0)
    p = Series60Parameters(
        length=L, beam=B, draft=T, cb=cb, nx=grid[0], nz=grid[1],
        transom_ratio=transom, cwp=cwp, lcb=lcb,
        bulb_length_ratio=bulb_length, bulb_area_ratio=bulb_area,
    )
    return generate_series60(p)


def _errors(model) -> dict[str, float | None]:
    return series60_analytic_relative_errors(model.parameters, model.hydrostatics)


def _finite(errors: dict[str, float | None]) -> dict[str, float]:
    return {k: v for k, v in errors.items() if isinstance(v, float) and math.isfinite(v)}


class TestScreenSweepOverTheDeclaredBox:
    """Every parent x transom in the box produces a comparison table with no holes."""

    @pytest.mark.parametrize("cb", PARENT_CB)
    @pytest.mark.parametrize("transom", TRANSOMS)
    def test_table_is_finite_and_within_the_screen_bar(self, cb: float, transom: float) -> None:
        errors = _errors(_model(cb, transom=transom))
        # ``nan`` is the failure mode that hid here: it is not a number, it is a
        # comparison that never happened.  ``None`` is the documented "no closed form".
        assert not any(isinstance(v, float) and math.isnan(v) for v in errors.values())
        compared = _finite(errors)
        for name in ALWAYS_COMPARED:
            assert name in compared, f"{name} was not compared at Cb={cb} transom={transom}"
        worst = max(compared, key=lambda k: abs(compared[k]))
        assert abs(compared[worst]) <= SCREEN_TOLERANCE, (
            f"Cb={cb} transom={transom}: {worst} is off by {compared[worst]:.4%}, "
            f"above the {SCREEN_TOLERANCE:.1%} screen"
        )

    @pytest.mark.parametrize("cb", PARENT_CB)
    @pytest.mark.parametrize("transom", TRANSOMS)
    def test_vcb_error_is_coupled_to_the_volume_error(self, cb: float, transom: float) -> None:
        compared = _finite(_errors(_model(cb, transom=transom)))
        volume = abs(compared["displacement_volume"])
        vcb = abs(compared["vertical_centre_buoyancy"])
        assert vcb <= volume + VCB_VOLUME_COUPLING, (
            f"Cb={cb} transom={transom}: VCB error {vcb:.4%} is not of the same order as "
            f"the volume error {volume:.4%} -- the two are not in the same frame"
        )


class TestConvergenceOnTheFiveParents:
    """The coarse screen is not a convergence claim; this is."""

    @pytest.mark.parametrize("cb", PARENT_CB)
    def test_fine_grid_is_inside_half_the_module_tolerance(self, cb: float) -> None:
        compared = _finite(_errors(_model(cb, grid=FINE_GRID)))
        worst = max(compared, key=lambda k: abs(compared[k]))
        assert abs(compared[worst]) <= FINE_TOLERANCE, (
            f"Cb={cb} at {FINE_GRID}: {worst} off by {compared[worst]:.4%}"
        )

    @pytest.mark.parametrize("cb", PARENT_CB)
    def test_refinement_reduces_the_error(self, cb: float) -> None:
        coarse = abs(_finite(_errors(_model(cb)))["displacement_volume"])
        fine = abs(_finite(_errors(_model(cb, grid=FINE_GRID)))["displacement_volume"])
        assert fine < coarse, (
            f"Cb={cb}: refining {SCREEN_GRID} -> {FINE_GRID} did not reduce the volume "
            f"error ({coarse:.4%} -> {fine:.4%}); the grid, not the surface, is wrong"
        )


class TestPublishedTargetsAcrossTheLadder:
    """``(cwp, lcb)`` targets are generator closures; the mesh must reproduce them."""

    @pytest.mark.parametrize("cb", PARENT_CB)
    def test_mesh_reproduces_the_targeted_waterplane_and_lcb(self, cb: float) -> None:
        model = _model(cb, cwp=PUBLISHED_CWP[cb], lcb=PUBLISHED_LCB)
        compared = _finite(_errors(model))
        assert abs(compared["waterplane_coefficient"]) <= CLOSURE_TOLERANCE
        lcb_error = abs(model.hydrostatics.longitudinal_centre_buoyancy / (PUBLISHED_LCB * L) - 1.0)
        assert lcb_error <= CLOSURE_TOLERANCE, (
            f"Cb={cb}: mesh LCB/L is {model.hydrostatics.longitudinal_centre_buoyancy / L:.5f}, "
            f"target {PUBLISHED_LCB}"
        )
        # The block coefficient is the volume target and must not move when the
        # waterplane is re-solved: the family is a Cb ladder first of all.
        assert abs(compared["block_coefficient"]) <= SCREEN_TOLERANCE

    @pytest.mark.parametrize("cb", PARENT_CB)
    def test_vcb_ratio_rises_with_block_coefficient(self, cb: float) -> None:
        """Fuller sections put the centre of buoyancy higher: physics, not algebra."""
        model = _model(cb)
        ratio = model.hydrostatics.vertical_centre_buoyancy / T + 1.0
        reference = analytic_series60_hydrostatics(model.parameters)["vertical_centre_buoyancy"] / T
        assert ratio == pytest.approx(reference, rel=0.03)
        assert 0.5 < ratio < 0.62


class TestSingleSourceOfTruthAndTheBulb:
    """The entry and the library must report the same table, and the bulb must not
    silently enter a hull-only comparison."""

    def test_cli_document_matches_the_library(self) -> None:
        from tensorcad.series60_cli import relative_errors

        model = _model(0.70, transom=0.2, grid=FINE_GRID)
        library = _errors(model)
        entry = relative_errors(model)
        assert set(entry) == set(library)
        for key, value in library.items():
            if value is None:
                assert entry[key] is None, f"{key}: entry fabricated {entry[key]!r}"
            else:
                assert entry[key] == pytest.approx(value, rel=1e-15, abs=1e-15)

    def test_bulb_leaves_the_hull_comparison_untouched(self) -> None:
        plain = _errors(_model(0.70, grid=SCREEN_GRID))
        with_bulb = _errors(_model(0.70, grid=SCREEN_GRID, bulb=(0.06, 0.05)))
        for name in ALWAYS_COMPARED:
            assert plain[name] == pytest.approx(with_bulb[name], abs=1e-12), (
                f"{name} moved when a bulb was added: the hull comparison picked up an "
                "appendage it has no mesh for"
            )
        assert not any(isinstance(v, float) and math.isnan(v) for v in with_bulb.values())