"""Topology, hydrostatics, continuity and defect-library tests for the pipeline.

Every test either proves a positive (the detector fires, the integral is exact) or a
*negative control* (a sub-tolerance gap is not a crack, a pristine hull raises no code,
a damaged hull cannot pass the gate).  The negative controls are the reason to believe
the positives.
"""
from __future__ import annotations

import pytest

from tensorcad.geomquality import baselines as baselines_module
from tensorcad.geomquality import continuity, defects
from tensorcad.geomquality.brep import audit_topology, orient_shell, sew
from tensorcad.geomquality.hydro import (
    convergence_study,
    facet_hydrostatics,
    patch_hydrostatics,
    richardson_gci,
)
from tensorcad.geomquality.pipeline import (
    CRITERIA,
    GateError,
    enforce,
    gate_single,
    verify_cad,
)


# --------------------------------------------------------------------------- #
# hydrostatics
# --------------------------------------------------------------------------- #
def test_wigley_patch_quadrature_is_exact() -> None:
    """Wigley III is a polynomial, so the patch route has no quadrature error left."""
    specimen = baselines_module.wigley_iii()
    hydro = patch_hydrostatics(
        specimen.brep, length=specimen.length, beam=specimen.beam, draft=specimen.draft
    )
    closed = specimen.reference["closed_form"]
    assert hydro.volume == pytest.approx(closed["volume"], rel=1e-12)
    assert hydro.waterplane_area == pytest.approx(closed["waterplane_area"], rel=1e-12)
    assert hydro.vcb_z == pytest.approx(closed["vcb_z"], abs=1e-12)
    assert hydro.lcb_x == pytest.approx(0.0, abs=1e-12)
    assert hydro.block_coefficient == pytest.approx(4.0 / 9.0, rel=1e-12)


def test_series60_closed_forms() -> None:
    """Cp, Cm, Cb and VCB/T are analytic for this parameterisation; check all four."""
    specimen = baselines_module.series60(cb=0.60)
    hydro = patch_hydrostatics(
        specimen.brep, length=specimen.length, beam=specimen.beam, draft=specimen.draft
    )
    closed = specimen.reference["closed_form"]
    assert hydro.volume == pytest.approx(closed["volume"], rel=1e-6)
    assert hydro.waterplane_area == pytest.approx(closed["waterplane_area"], rel=1e-6)
    assert hydro.vcb_z == pytest.approx(closed["vcb_z"], rel=1e-6)
    assert hydro.block_coefficient == pytest.approx(0.60, rel=1e-6)


def test_kvlcc2_class_reproduces_published_particulars() -> None:
    specimen = baselines_module.kvlcc2_class()
    hydro = patch_hydrostatics(
        specimen.brep, length=specimen.length, beam=specimen.beam, draft=specimen.draft
    )
    assert hydro.block_coefficient == pytest.approx(0.8098, rel=1e-4)
    assert abs(hydro.lcb_x) / specimen.length == pytest.approx(0.0348, abs=1e-4)
    assert hydro.lcb_x < 0  # aft of midship, as published
    assert hydro.volume > 300_000


def test_negative_control_beam_scaled_is_caught_by_cb() -> None:
    """A 5 % narrow hull integrates to the wrong Cb: that must not pass as clean."""
    specimen = baselines_module.wigley_iii()
    declared = {
        "length": specimen.length,
        "beam": specimen.beam,
        "draft": specimen.draft,
        "block_coefficient": 4.0 / 9.0,
    }
    damaged = defects.scale_brep(specimen.brep, 0.95, axes=(1,), source="narrow")
    detection = defects.detect_all(damaged, declared=declared)
    assert "CB_MISMATCH" in detection["codes"]
    assert detection["watertight"] is True  # topology is fine -- the ship is wrong
    clean = defects.detect_all(specimen.brep, declared=declared)
    assert clean["codes"] == []


def test_facet_route_error_decreases_and_converges_second_order() -> None:
    specimen = baselines_module.wigley_iii()
    reference = specimen.reference["closed_form"]["volume"]
    errors = []
    for u_samples, v_samples in ((12, 6), (24, 12), (48, 24)):
        table = facet_hydrostatics(specimen.brep, u_samples=u_samples, v_samples=v_samples)
        errors.append(abs(table.volume - reference) / reference)
    assert errors[0] > errors[1] > errors[2], errors
    study = convergence_study(
        specimen.brep, levels=((12, 6), (24, 12), (48, 24)), reference=reference
    )
    assert study["observed_order"] > 1.9
    assert study["gci_percent"] < 1.0


def test_physical_refinement_restores_second_order_on_a_graded_patch() -> None:
    """Parameter-space refinement is not physical refinement.

    The Series 60 patch maps ``z = -T + T*v**24``, so uniform ``v`` samples put almost
    every node near the keel; refining that sequence measures the wrong thing and the
    observed order collapses.  Equal-arc-length sampling refines the surface itself.
    """
    specimen = baselines_module.series60(cb=0.60)
    reference = specimen.reference["volume"]
    levels = ((8, 4), (16, 8), (32, 16))
    physical = convergence_study(
        specimen.brep, levels=levels, reference=reference, physical=True
    )
    parameter = convergence_study(
        specimen.brep, levels=levels, reference=reference, physical=False
    )
    assert physical["observed_order"] > 1.9
    assert physical["observed_order"] > parameter["observed_order"]
    assert physical["gci_percent"] < parameter["gci_percent"]


def test_arc_length_parameterisation_is_monotone_and_end_pointed() -> None:
    from tensorcad.geomquality.hydro import arc_length_parameters

    specimen = baselines_module.series60(cb=0.60)
    hull = [face for face in specimen.brep.faces if face.name.endswith("hull_half")][0]
    parameters = arc_length_parameters(hull.surface, 1, 16, probes=65)
    assert parameters[0] == pytest.approx(0.0, abs=1e-12)
    assert parameters[-1] == pytest.approx(1.0, abs=1e-12)
    assert all(b > a for a, b in zip(parameters, parameters[1:]))
    # Most of the physical surface is packed into the last part of the parameter range,
    # so equal-arc sampling puts the median node at v > 0.8, not at v = 0.5.
    median_depth = -specimen.draft + specimen.draft * parameters[8] ** 24
    assert parameters[8] > 0.8
    assert median_depth < -0.75 * specimen.draft


def test_richardson_reports_exact_discretisation_instead_of_fitting() -> None:
    report = richardson_gci([2.0, 2.0, 2.0])
    assert report["converged"] is True
    assert report["observed_order"] is None
    assert report["gci_percent"] == 0.0


# --------------------------------------------------------------------------- #
# topology
# --------------------------------------------------------------------------- #
def test_clean_hulls_are_watertight_with_score_one() -> None:
    for name, factory in baselines_module.specimens().items():
        specimen = factory()
        audit = audit_topology(specimen.brep)
        assert audit.watertight is True, name
        assert audit.score == 1.0
        assert audit.defects == ()
        assert audit.signed_volume > 0


def test_sub_tolerance_gap_is_not_a_crack() -> None:
    """A 1e-7 m mismatch is weld noise; a 5e-5 m mismatch (5e-6 L) is a crack.

    This is the negative control for the crack detector: without it, the detector would
    be indistinguishable from one that reports every pair of nearby curves.
    """
    specimen = baselines_module.wigley_iii()
    subtle = defects._seam_gap_offset(specimen.brep, delta=1e-7)
    visible = defects._seam_gap_offset(specimen.brep, delta=5e-5)
    subtle_codes = {defect.code for defect in audit_topology(subtle).defects}
    visible_audit = audit_topology(visible)
    visible_codes = {defect.code for defect in visible_audit.defects}
    assert "SEAM_GAP" not in subtle_codes
    assert "SEAM_GAP" in visible_codes
    assert visible_audit.watertight is False
    assert visible_audit.max_gap > 1e-5


def test_degenerate_edges_are_legal_seams_not_micro_edges() -> None:
    """The bow/stern of a hull collapse to a point; that is not a tolerance violation."""
    spec = baselines_module.series60(cb=0.60)
    audit = audit_topology(spec.brep)
    assert audit.watertight is True
    assert "MICRO_EDGE" not in {defect.code for defect in audit.defects}
    assert audit.degenerate_edges  # recorded, just not fatal


def test_orientation_propagation_repairs_a_randomly_oriented_shell() -> None:
    """IGES gives no orientation; the kernel has to derive one and pick the outward side."""
    specimen = baselines_module.wigley_iii()
    random_shell = sew(
        [face.surface for face in specimen.brep.faces],
        names=[face.name for face in specimen.brep.faces],
        same_sense=[True] * len(specimen.brep.faces),
        wetted=[face.wetted for face in specimen.brep.faces],
        source="randomly oriented",
    )
    assert audit_topology(random_shell).watertight is False
    repaired = orient_shell(random_shell)
    audit = audit_topology(repaired)
    assert audit.watertight is True
    assert audit.signed_volume > 0


# --------------------------------------------------------------------------- #
# continuity
# --------------------------------------------------------------------------- #
def test_continuity_thresholds_calibrate_and_score_above_nine_tenths() -> None:
    evidence = continuity.continuity_evidence(cases=90)
    assert evidence["calibration_cases"] >= 30
    assert evidence["test_cases"] >= 20
    assert evidence["test"]["overall"]["f1"] > 0.9
    assert evidence["test"]["overall"]["precision"] > 0.9
    # The thresholds must actually separate the classes, not just fit the data.
    separation = evidence["separation"]
    assert separation["curvature_jump_min_G1"] > separation["curvature_jump_max_G2"]
    assert separation["normal_angle_min_G0_deg"] > 1.0


def test_g1_seam_is_not_g2_and_g2_seam_is_not_flagged() -> None:
    for kind, expected in (("G2", "G2"), ("G1", "G1"), ("G0", "G0")):
        import random

        first, second, _ = continuity.make_seam_case(kind, random.Random(11))
        from tensorcad.geomquality.brep import Face

        metrics = continuity.seam_metrics(
            Face(index=0, surface=first, boundary=()),
            Face(index=1, surface=second, boundary=()),
            samples=9,
            scale=1.0,
        )
        assert (
            continuity.classify(
                metrics,
                angle_limit_deg=0.05,
                curvature_limit=0.01,
            )
            == expected
        ), (kind, metrics)


def test_designed_creases_are_not_continuity_defects() -> None:
    """The keel knife edge and the deck chine are G0 by design; they must read KNUCKLE."""
    specimen = baselines_module.wigley_iii()
    seams = continuity.audit_seams(specimen.brep)
    assert seams
    assert {seam.classification for seam in seams} <= {"KNUCKLE", "G2"}
    assert not [seam for seam in seams if seam.classification in {"G0", "G1"}]


# --------------------------------------------------------------------------- #
# defect library and gate
# --------------------------------------------------------------------------- #
def test_defect_library_detects_every_case_without_false_alarms() -> None:
    """The library runs on the fast reference hull; the full evidence bundle covers all three.

    (A Series 60 case costs seconds because its vertical span is degree 24 -- exact
    quadrature to 3e-8, which is the point of the pipeline but not of a unit test.)
    """
    name = "wigley_iii"
    specimen = baselines_module.specimens()[name]()
    declared = {
        "length": specimen.length,
        "beam": specimen.beam,
        "draft": specimen.draft,
        "block_coefficient": specimen.particulars["block_coefficient"],
    }
    result = defects.run_library(specimen.brep, declared=declared)
    assert result["cases"] >= 10
    assert result["detected"] == result["cases"], [
        row for row in result["rows"] if not row["detected"]
    ]
    assert result["false_alarm_codes"] == []
    assert result["baseline_watertight"] is True
    for row in result["rows"]:
        assert row["repair"], row["case"]  # every code carries a repair hint


def test_gate_accepts_a_clean_hull_and_rejects_a_damaged_one() -> None:
    specimen = baselines_module.wigley_iii()
    declared = {
        "length": specimen.length,
        "beam": specimen.beam,
        "draft": specimen.draft,
        "block_coefficient": 4.0 / 9.0,
    }
    clean = verify_cad(specimen.brep, declared=declared, reference=specimen.reference)
    enforce(gate_single(clean))

    holed = defects.DEFECT_CASES[0].inject(specimen.brep)
    damaged = verify_cad(holed, declared=declared, reference=specimen.reference)
    assert damaged["watertight"] is False
    assert damaged["codes"]
    with pytest.raises(GateError):
        enforce(gate_single(damaged))


def test_verify_cad_reports_timings_and_criteria_numbers() -> None:
    specimen = baselines_module.wigley_iii()
    declared = {"length": 10.0, "beam": 1.0, "draft": 0.625, "block_coefficient": 4.0 / 9.0}
    report = verify_cad(
        specimen.brep,
        declared=declared,
        reference=specimen.reference,
        convergence_levels=((8, 4), (16, 8), (32, 16)),
    )
    assert report["seconds"] < 10.0
    assert report["comparison"]["volume"]["patch_relative_error"] < 1e-10
    assert report["convergence"]["observed_order"] > 1.9
    assert report["brep"]["faces"] == 4
    assert set(CRITERIA) >= {"C1_wigley_volume_error", "G1_watertight"}
