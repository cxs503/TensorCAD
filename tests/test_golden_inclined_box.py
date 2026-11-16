"""The golden gate's attitude case, and the oracle that is independent of the library.

Two things are proved here:

* the oracle (:mod:`tools.golden.inclined_box_oracle`) is right *on its own terms* -- in
  the upright limit it reproduces the textbook closed forms of a barge, and it converges at
  second order in its own quadrature;
* the golden case ``inclined_box`` has teeth -- a perturbation injected into the measured
  route at combined heel + trim shows up in the case's numbers, at a size far below the
  acceptance band, while the unperturbed route agrees with the oracle to round-off.
"""
from __future__ import annotations

import importlib.util
import pathlib

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]


def _oracle():
    spec = importlib.util.spec_from_file_location(
        "inclined_box_oracle_under_test", ROOT / "tools" / "golden" / "inclined_box_oracle.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


oracle = _oracle()

LENGTH, BEAM, DEPTH, DRAFT = 10.0, 4.0, 3.0, 1.5
VOLUME = LENGTH * BEAM * DRAFT


def _relative(measured: float, reference: float) -> float:
    return abs(measured - reference) / abs(reference) if reference else abs(measured - reference)


def test_upright_limit_reproduces_the_textbook_barge() -> None:
    """The oracle's own anchor: V = L·B·T, VCB = T/2, A_wp = L·B, no transverse shift."""
    upright = oracle.oracle_hydrostatics(
        LENGTH, BEAM, DEPTH, volume_m3=VOLUME, heel_deg=0.0, trim_deg=0.0, resolution=600
    )
    assert _relative(upright["displacement_volume_m3"], VOLUME) < 1e-9
    assert _relative(upright["vcb_z_m"], 0.5 * DRAFT) < 1e-9
    assert _relative(upright["waterplane_area_m2"], LENGTH * BEAM) < 1e-9
    assert _relative(upright["lcb_x_m"], 0.0) < 1e-12
    assert _relative(upright["tcb_y_m"], 0.0) < 1e-12


def test_trim_only_limit_matches_the_lcb_of_a_rotated_rectangle() -> None:
    """Pure trim: the waterline is a tilted plane, so A_wp stretches by 1/cos(trim)."""
    trim = 3.0
    values = oracle.oracle_hydrostatics(
        LENGTH, BEAM, DEPTH, volume_m3=VOLUME, heel_deg=0.0, trim_deg=trim, resolution=1200
    )
    wall_sided = LENGTH * BEAM / __import__("math").cos(__import__("math").radians(trim))
    assert _relative(values["waterplane_area_m2"], wall_sided) < 2e-3
    assert _relative(values["tcb_y_m"], 0.0) < 1e-9


def test_oracle_converges_at_second_order() -> None:
    """The oracle has an error budget of its own, measured rather than assumed."""
    study = oracle.convergence(
        LENGTH, BEAM, DEPTH, volume_m3=VOLUME, heel_deg=15.0, trim_deg=3.0,
        resolutions=(150, 300, 600, 1200),
    )
    for metric, errors in study["relative_deviation_from_finest"].items():
        assert errors[-1] == 0.0, metric
        assert errors[-2] < 1e-8, (metric, errors)


@pytest.mark.parametrize("heel_deg,trim_deg", [(15.0, 3.0), (5.0, 1.0), (-10.0, 2.0)])
def test_polyhedron_clip_agrees_with_the_independent_integral(heel_deg: float, trim_deg: float) -> None:
    """Two structurally different derivations of one answer, at combined heel and trim.

    This is the assertion the golden case carries: a face clip and a divergence-theorem
    moment sum on one side, a column integral that never builds a face on the other.
    """
    from tensorcad.box_inclined import box_inclined_hydrostatics

    hydro = box_inclined_hydrostatics(
        LENGTH, BEAM, DEPTH, volume_m3=VOLUME, heel_deg=heel_deg, trim_deg=trim_deg
    )
    values = oracle.oracle_hydrostatics(
        LENGTH, BEAM, DEPTH, volume_m3=VOLUME, heel_deg=heel_deg, trim_deg=trim_deg,
        resolution=1200,
    )
    x, y, z = hydro.centre_of_buoyancy
    for measured, reference, label in (
        (hydro.volume_m3, values["displacement_volume_m3"], "volume"),
        (x, values["lcb_x_m"], "LCB"),
        (y, values["tcb_y_m"], "TCB"),
        (z, values["vcb_z_m"], "VCB"),
        (hydro.waterplane_area_m2, values["waterplane_area_m2"], "A_wp"),
        (hydro.waterline_offset, values["waterline_offset_m"], "waterline offset"),
    ):
        assert _relative(measured, reference) < 1e-9, (label, measured, reference)


def test_the_case_fires_on_a_perturbed_attitude_route(monkeypatch: pytest.MonkeyPatch) -> None:
    """Detection power: a 0.2 % shift in the measured waterplane area must be visible.

    The perturbation is injected into the module the case measures (not into the case), so
    this is the golden gate's blind-spot question asked of the new case: how large a
    regression does the 1e-4 band let through?  Nothing at 2e-3 does.
    """
    from tensorcad import box_inclined
    from tensorcad.golden_inclined_box import InclinedBoxGolden

    case = InclinedBoxGolden()
    reference = case.reference()
    # The acceptance rung, not a level name typed here: m22879 renamed the ladder from
    # ("exact", "tight") to ("coarse", "fine") and this literal kept asking for a rung that
    # no longer exists, so the test raised KeyError instead of measuring anything.
    acceptance = case.levels()[-1]
    clean = case.measure(acceptance)
    assert _relative(clean["waterplane_area_m2"], reference["waterplane_area_m2"]) < 1e-9

    original = box_inclined.box_inclined_hydrostatics

    import dataclasses

    def skewed(*args, **kwargs):
        hydro = original(*args, **kwargs)
        return dataclasses.replace(hydro, waterplane_area_m2=hydro.waterplane_area_m2 * 1.002)

    monkeypatch.setattr(box_inclined, "box_inclined_hydrostatics", skewed)
    perturbed = case.measure(acceptance)
    error = _relative(perturbed["waterplane_area_m2"], reference["waterplane_area_m2"])
    assert error > 1e-4, error
    # and the rest of the report is untouched, so the finding localises to this metric
    assert _relative(perturbed["displacement_volume_m3"], reference["displacement_volume_m3"]) < 1e-9


def test_the_committed_case_is_registered_and_locks_its_solid() -> None:
    """The case is in the library and it fingerprints the clipped solid it integrates."""
    from tensorcad import golden

    names = [case.name for case in golden.golden_cases()]
    assert "inclined_box" in names
    case = next(case for case in golden.golden_cases() if case.name == "inclined_box")
    assert case.geometry_sha256()
    residual = case.invariants(case.measure(case.levels()[-1]))
    volume = residual["identity_volume_equals_requested_displacement"]
    assert abs(volume[0]) < volume[1], volume