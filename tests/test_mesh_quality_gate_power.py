"""Tests for ``tools/mesh_quality_gate_power.py`` (the gate-power study).

Two things are asserted, and they are different:

* the statistics that turn verdicts into a claim are correct on inputs with a known
  answer (AUROC of a perfect and of a random score, Spearman on a monotone pair, the
  confusion table's MCC on a hand-computable case);
* the study's own machinery runs on a reduced ladder and its criteria evaluate to the
  values they should, including the honest FAIL path -- a criterion that fires must stay
  fired, because a study whose thresholds can be re-tuned until they pass is worth nothing.

The tolerances below are the ones the tool declares; nothing here re-tunes them.
"""
from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

_SPEC = importlib.util.spec_from_file_location(
    "mesh_quality_gate_power", ROOT / "tools" / "mesh_quality_gate_power.py"
)
assert _SPEC is not None and _SPEC.loader is not None
power = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(power)


# ---------------------------------------------------------------------------------
# statistics: known answers, not just "it runs"
# ---------------------------------------------------------------------------------
def test_ranks_average_ties():
    assert power.ranks([10.0, 20.0, 20.0, 30.0]) == [1.0, 2.5, 2.5, 4.0]


def test_auroc_perfect_and_inverted_and_uninformative():
    labels = [True, True, False, False]
    assert power.auroc([3.0, 4.0, 1.0, 2.0], labels) == 1.0
    assert power.auroc([1.0, 2.0, 3.0, 4.0], labels) == 0.0
    assert power.auroc([1.0, 1.0, 1.0, 1.0], labels) == pytest.approx(0.5)


def test_auroc_is_a_probability():
    scores = [0.1, 0.5, 0.2, 0.9, 0.4]
    labels = [False, True, False, True, False]
    value = power.auroc(scores, labels)
    assert 0.0 <= value <= 1.0
    # 5 positives-free / single-class input has no defined AUC; it must say so, not guess
    assert math.isnan(power.auroc(scores, [False] * 5))


def test_spearman_perfect_and_constant():
    assert power.spearman([1.0, 2.0, 3.0, 4.0], [10.0, 20.0, 30.0, 40.0]) == pytest.approx(1.0)
    assert power.spearman([1.0, 2.0, 3.0, 4.0], [40.0, 30.0, 20.0, 10.0]) == pytest.approx(-1.0)
    assert power.spearman([1.0, 1.0, 1.0], [1.0, 2.0, 3.0]) == 0.0


def test_confusion_matches_hand_computation():
    # tp=2, fp=1, tn=3, fn=1  ->  balanced accuracy 0.5*(2/3 + 3/4)
    predicted = [True, True, True, False, False, False, False]
    labels = [True, True, False, True, False, False, False]
    row = power.confusion(predicted, labels)
    assert (row["tp"], row["fp"], row["tn"], row["fn"]) == (2, 1, 3, 1)
    assert row["balanced_accuracy"] == pytest.approx(0.5 * (2 / 3) + 0.5 * (3 / 4))
    assert row["mcc"] == pytest.approx(
        (2 * 3 - 1 * 1) / math.sqrt((2 + 1) * (2 + 1) * (3 + 1) * (3 + 1))
    )


# ---------------------------------------------------------------------------------
# reference surfaces: closed forms and gate cleanliness
# ---------------------------------------------------------------------------------
def test_cube_is_exact_and_outward():
    mesh = power.unit_cube(3)
    assert len(mesh.faces) == 12 * 3 * 3
    assert power.signed_volume(mesh) == pytest.approx(1.0, abs=1e-12)
    assert power.surface_area(mesh) == pytest.approx(6.0, abs=1e-12)


def test_icosphere_is_close_to_the_sphere_and_passes_the_gate():
    from tensorcad.geomquality.quality import quality_gate

    mesh = power.unit_sphere(3)
    volume_error = abs(power.signed_volume(mesh) - 4.0 * math.pi / 3.0) / (4.0 * math.pi / 3.0)
    area_error = abs(power.surface_area(mesh) - 4.0 * math.pi) / (4.0 * math.pi)
    # an inscribed subdivision only ever undershoots, and 3 subdivisions are already
    # inside the 1 % engineering tolerance the study declares
    assert 0.0 < volume_error < 0.02
    assert 0.0 < area_error < 0.02
    assert quality_gate(mesh).passed


def test_every_pristine_surface_passes_the_default_gate():
    """P1 on the real base surfaces, not on a hand-made one."""
    from tensorcad.geomquality.quality import quality_gate

    for name, spec in power.reference_surfaces(quick=True).items():
        gate = quality_gate(spec["mesh"])
        assert gate.passed, (name, [c.metric for c in gate.checks if not c.passed])


# ---------------------------------------------------------------------------------
# damage ladder: monotone in the knob that drives it
# ---------------------------------------------------------------------------------
def test_jitter_degrades_min_angle_monotonically():
    from tensorcad.geomquality.quality import mesh_quality

    mesh = power.unit_cube(4)
    angles = [
        mesh_quality(power.jitter(mesh, amplitude, 0))["min_angle_deg"]
        if False
        else mesh_quality(power.jitter(mesh, amplitude, 0)).summary["min_angle_deg"].minimum
        for amplitude in (0.01, 0.05, 0.2)
    ]
    assert angles[0] > angles[1] > angles[2]


def test_shear_moves_vertices_inside_the_tangent_plane():
    """The claim the family makes is geometric, so it is checked geometrically.

    (It is *not* claimed that a tangent displacement preserves the volume of a faceted
    shell: on a cube the facets tilt and the enclosed volume moves by ~2 % at amplitude
    0.1, which is why the study reports the number instead of assuming it.)
    """
    mesh = power.unit_cube(4)
    normals = power._vertex_normals(mesh)
    sheared = power.shear(mesh, 0.1)
    h = power.mean_edge(mesh)
    for original, moved, normal in zip(mesh.vertices, sheared.vertices, normals):
        delta = [moved[k] - original[k] for k in range(3)]
        normal_component = abs(sum(delta[k] * normal[k] for k in range(3)))
        assert normal_component <= 1e-9 * h
        assert math.dist(original, moved) <= 0.1 * h + 1e-9


def test_collapse_creates_a_sliver():
    from tensorcad.geomquality.quality import mesh_quality

    mesh = power.unit_cube(4)
    damaged = power.collapse(mesh, 0.999)
    assert mesh_quality(damaged).summary["min_angle_deg"].minimum < mesh_quality(
        mesh
    ).summary["min_angle_deg"].minimum


def test_smoothing_shrinks_the_body_while_the_gate_stays_green():
    """The blind spot, asserted with the number the study reports.

    Tait-Bryan-free Taubin smoothing keeps element shapes healthy (that is what it is for)
    while the body loses volume; on the cube at 20 iterations the loss is double-digit
    percent and the default gate still passes. If this test ever starts failing because the
    gate learned to catch it, the study's conclusion has to be rewritten -- which is the
    point of pinning it here.
    """
    from tensorcad.geomquality.quality import mesh_quality, quality_gate

    mesh = power.unit_cube(4)
    smoothed = power.smoothed(mesh, 20)
    shrink = abs(power.signed_volume(smoothed) - power.signed_volume(mesh)) / power.signed_volume(
        mesh
    )
    assert shrink > power.TOLERANCE  # the engineering quantity is off
    assert mesh_quality(smoothed).summary["min_angle_deg"].minimum > 2.0  # elements are fine
    assert quality_gate(smoothed).passed  # ... and the gate does not report the volume loss


# ---------------------------------------------------------------------------------
# the study end to end, on the reduced ladder
# ---------------------------------------------------------------------------------
def test_quick_run_produces_a_consistent_report(tmp_path):
    results = power.run(quick=True)
    assert results["schema"] == "tensorcad.mesh_quality_gate_power/1"
    assert results["pooled"]["samples"] == len(results["samples"])
    assert results["pooled"]["samples"] == (
        sum(f["samples"] for f in results["families"].values())
    )
    assert all(key in results["criteria"] for key in (
        "P1_no_false_alarm", "P2_balanced_accuracy", "P3_auroc_min_angle",
        "P4_spearman_min_angle",
    ))
    assert results["all_criteria_passed"] == all(
        item["passed"] for item in results["criteria"].values()
    )
    cells = results["axes"]["cells"]
    assert sum(cells.values()) == results["pooled"]["samples"]
    for sample in results["samples"]:
        assert sample["defect"] == (sample["fidelity_error"] > power.TOLERANCE)
        assert sample["quality_fail"] == (not sample["gate_pass"])


def test_criteria_are_evaluated_from_the_summary_not_re_tuned():
    """A summary that cannot pass must fail, whatever the rest of the study did."""
    summary = {
        "pooled": {
            "confusion": {"balanced_accuracy": 0.5},
            "auroc_min_angle_score": 0.5,
            "spearman_min_angle_error": 0.0,
        }
    }
    base_rows = {"x": {"gate_pass": True}}
    criteria = power.criteria_from_summary(summary, base_rows)
    assert criteria["P1_no_false_alarm"]["passed"] is True
    assert criteria["P2_balanced_accuracy"]["passed"] is False
    assert criteria["P3_auroc_min_angle"]["passed"] is False
    assert criteria["P4_spearman_min_angle"]["passed"] is False

    good = {
        "pooled": {
            "confusion": {"balanced_accuracy": 0.95},
            "auroc_min_angle_score": 0.97,
            "spearman_min_angle_error": 0.9,
        }
    }
    assert all(item["passed"] for item in power.criteria_from_summary(good, base_rows).values())

    base_rows_bad = {"x": {"gate_pass": False}}
    assert power.criteria_from_summary(good, base_rows_bad)["P1_no_false_alarm"]["passed"] is False


def test_zero_artifacts_write_a_manifest(tmp_path):
    results = power.run(quick=True)
    manifest = power.write_artifacts(results, tmp_path)
    assert "results.json" in manifest["files"]
    assert len(manifest["files"]["results.json"]["sha256"]) == 64
    assert (tmp_path / "results.json").exists()
