"""Out-of-distribution audit of the G0/G1/G2 seam classifier (``continuity_ood``).

The production evidence for the classifier is a held-out split of *one* generator
(planar cubic plates, straight seam). These tests pin the two things that split cannot
show: the regime in which the "relative" curvature criterion is really absolute, and the
resolution at which a localized defect disappears between the 17 sample points -- plus
the adaptive sampler that closes the second one without disturbing the first.
"""
from __future__ import annotations

import random

import pytest

from tensorcad.geomquality import baselines
from tensorcad.geomquality import continuity
from tensorcad.geomquality import continuity_ood as ood
from tensorcad.geomquality.brep import Face

#: Dense enough that the 1/64-wide defect is crossed by dozens of samples, cheap enough
#: for a unit test (the full study in ``examples/geomquality/audit_continuity_ood.py``
#: uses 2049).
DENSE = 1025


def _metrics(face_a: Face, face_b: Face, samples: int, *, scale: float = 1.0) -> dict:
    return continuity.seam_metrics(face_a, face_b, samples=samples, scale=scale)


# --------------------------------------------------------------------------- #
# the constructors really are what the study assumes
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("curvature", [0.05, 1.0, 7.5])
def test_curved_plate_pair_has_the_curvature_jump_it_was_built_with(curvature: float) -> None:
    """``curved_plate`` + ``straight_plate`` is G1 with |dkappa| == curvature, exactly."""
    flat = ood._face(ood.straight_plate(1.0, 1.0, incoming=True), 0)
    bent = ood._face(ood.curved_plate(curvature, 1.0, 1.0, incoming=False), 1)
    metrics = _metrics(flat, bent, 65)
    assert metrics["max_normal_angle_deg"] == pytest.approx(0.0, abs=1e-9)
    assert metrics["max_curvature_jump_absolute"] == pytest.approx(curvature, rel=1e-9)
    assert metrics["max_gap"] == pytest.approx(0.0, abs=1e-12)


def test_seam_metrics_reports_the_absolute_jump_on_both_branches() -> None:
    """``max_curvature_jump_absolute`` exists whether or not a seam was found."""
    flat = ood._face(ood.straight_plate(1.0, 1.0, incoming=True), 0)
    bent = ood._face(ood.curved_plate(1.0, 1.0, 1.0, incoming=False), 1)
    assert "max_curvature_jump_absolute" in _metrics(flat, bent, 9)
    # Two plates a metre apart share no isocurve: the no-seam branch, same keys.
    far = ood._face(ood.straight_plate(1.0, 1.0, incoming=False), 1)
    moved = Face(
        index=1,
        surface=continuity.NURBSSurface(
            degree=far.surface.degree,
            control_net=tuple(
                tuple((x + 5.0, y, z) for x, y, z in row) for row in far.surface.control_net
            ),
            knots=far.surface.knots,
        ),
        boundary=(),
    )
    assert "max_curvature_jump_absolute" in _metrics(flat, moved, 9)


# --------------------------------------------------------------------------- #
# study 1: relative above the floor, absolute below it
# --------------------------------------------------------------------------- #
def test_thirty_percent_mismatch_is_only_visible_above_the_scale_floor() -> None:
    """The docstring's "curvatures off by 30 %" needs kappa * scale > 0.1/0.3, not just 30 %."""
    rows = ood.curvature_scale_sweep(mismatch=0.3, levels=(0.1, 0.3, 1.0, 3.0))
    by_level = {round(row["level_times_scale"], 6): row for row in rows}
    # Same relative mismatch (30 %), same construction: only the absolute level changes.
    assert by_level[3.0]["max_curvature_jump"] == pytest.approx(0.3, rel=1e-6)
    assert by_level[3.0]["classification"] == "G1"
    assert by_level[0.1]["max_curvature_jump"] < continuity.DEFAULT_CURVATURE_LIMIT
    assert by_level[0.1]["classification"] == "G2"
    # ... and the absolute jump is 30 % of the level in both, i.e. the floor is the cause.
    assert by_level[0.1]["max_curvature_jump_absolute"] == pytest.approx(0.03, rel=1e-6)


def test_full_mismatch_switches_on_at_the_curvature_limit() -> None:
    rows = ood.curvature_scale_sweep(mismatch=1.0, levels=(1e-3, 0.03, 0.3, 3.0))
    verdicts = {round(row["level_times_scale"], 6): row["classification"] for row in rows}
    assert verdicts[1e-3] == "G2"  # a 100 % jump, invisible: below 1/scale
    assert verdicts[0.03] == "G2"
    assert verdicts[0.3] == "G1"
    assert verdicts[3.0] == "G1"
    # Below the floor the ratio is the absolute jump times the scale -- that is the model.
    low = next(row for row in rows if row["level_times_scale"] == 1e-3)
    assert low["max_curvature_jump"] == pytest.approx(1e-3, rel=1e-6)


# --------------------------------------------------------------------------- #
# study 2: a localized defect disappears between the sample points
# --------------------------------------------------------------------------- #
def _localized(width: float, centre: float = 0.53125, curvature: float = 2.0):
    face_a, face_b, info = ood.localized_curvature_defect(
        curvature, width, centre=centre
    )
    return face_a, face_b, info


def test_a_large_localised_defect_is_reported_clean_at_the_default_resolution() -> None:
    """kappa * scale = 2, i.e. 20x the limit, still reads G2 when it fits between samples."""
    face_a, face_b, info = _localized(1 / 64)
    truth = _metrics(face_a, face_b, DENSE)
    coarse = _metrics(face_a, face_b, ood.PRODUCTION_SAMPLES)
    assert truth["max_curvature_jump_absolute"] > 1.0  # the defect is real and severe
    assert continuity.classify(truth) == "G1"
    assert continuity.classify(coarse) == "G2"  # ... and the default grid misses it
    capture = coarse["max_curvature_jump"] / truth["max_curvature_jump"]
    assert capture < 0.1
    assert info["width"] == pytest.approx(1 / 64)


def test_the_same_defect_on_a_sample_point_is_found() -> None:
    """Control: centring the identical defect on a grid point restores the verdict.

    So the miss above is aliasing, not a defect too small to matter.
    """
    face_a, face_b, _ = _localized(1 / 64, centre=0.5)
    truth = _metrics(face_a, face_b, DENSE)
    coarse = _metrics(face_a, face_b, ood.PRODUCTION_SAMPLES)
    assert continuity.classify(truth) == "G1"
    assert continuity.classify(coarse) == "G1"
    assert coarse["max_curvature_jump"] / truth["max_curvature_jump"] > 0.9


def test_capture_ratio_collapses_as_the_defect_narrows() -> None:
    rows = ood.sampling_alias_sweep(
        widths=(1 / 8, 1 / 32, 1 / 64), centres=(0.53125,), dense=DENSE
    )
    captures = [row["coarse_capture"] for row in rows]
    assert captures == sorted(captures, reverse=True)
    assert captures[0] > 0.9
    assert captures[-1] < 0.1
    assert all(row["adaptive_recovered"] for row in rows)


# --------------------------------------------------------------------------- #
# the fix
# --------------------------------------------------------------------------- #
def test_adaptive_sampling_recovers_the_missed_defect() -> None:
    face_a, face_b, _ = _localized(1 / 64)
    adaptive = continuity.seam_metrics_adaptive(
        face_a, face_b, samples=ood.PRODUCTION_SAMPLES, max_samples=513, scale=1.0
    )
    assert adaptive["converged"] is True
    assert ood.PRODUCTION_SAMPLES < adaptive["samples_used"] <= 513
    assert continuity.classify(adaptive) == "G1"
    assert adaptive["max_curvature_jump"] > 0.5


def test_adaptive_sampling_leaves_the_calibrated_family_alone() -> None:
    """No verdict on the calibration family may move, or the limits would need refitting."""
    report = ood.adaptive_regression(cases=18, max_samples=513)
    assert report["cases"] == 18
    assert report["changed"] == 0
    # A clean seam settles at the second resolution: adaptivity is nearly free there.
    assert report["samples_used_max"] <= 64


def test_audit_seams_refine_agrees_with_the_fixed_grid_on_a_fair_hull() -> None:
    specimen = baselines.wigley_iii()
    fixed = continuity.audit_seams(specimen.brep)
    refined = continuity.audit_seams(specimen.brep, refine=True, max_samples=513)
    assert fixed and len(fixed) == len(refined)
    assert [seam.classification for seam in fixed] == [
        seam.classification for seam in refined
    ]
    assert all(seam.converged for seam in refined)
    assert all(seam.samples > ood.PRODUCTION_SAMPLES for seam in refined)


def test_adaptive_measurement_of_a_g0_seam_still_reports_g0() -> None:
    """Angle discontinuities are already resolved at 17 points; adaptivity must not blur them."""
    surface_a, surface_b, _ = continuity.make_seam_case("G0", random.Random(5))
    adaptive = continuity.seam_metrics_adaptive(
        Face(index=0, surface=surface_a, boundary=()),
        Face(index=1, surface=surface_b, boundary=()),
        samples=ood.PRODUCTION_SAMPLES,
        max_samples=513,
        scale=1.0,
    )
    assert continuity.classify(adaptive) == "G0"
    assert adaptive["max_normal_angle_deg"] > continuity.DEFAULT_ANGLE_LIMIT_DEG