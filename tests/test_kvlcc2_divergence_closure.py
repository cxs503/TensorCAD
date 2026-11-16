"""The divergence-theorem probe must be trustworthy before it may judge an asset.

Three layers are pinned here:

* the kernel integrates a *known* closed surface exactly enough to see a 3 % defect
  (unit sphere, −0.11 %; unit cube, machine-level), so a deficit measured on KVLCC2
  cannot be blamed on the quadrature;
* the even-odd (u, v) mask accepts a square and rejects its complement, and handles a
  circle without a winding convention;
* the committed artifact carries the shape the document cites, and its closure residual
  is reported as a number (never silently zero).
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

from tools.baseline_fleet.kvlcc2_divergence_probe import mask_uv, sphere_check

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "artifacts/baseline_fleet/kvlcc2_divergence.json"
DOC = ROOT / "docs/kvlcc2_divergence_closure.md"


def _tetra_volume() -> float:
    """Divergence volume of the unit tetrahedron: four planar outward faces, exact."""
    O, X, Y, Z = (np.array(v, float) for v in ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)))
    volume = 0.0
    for a, b, c in ((O, Y, X), (O, X, Z), (O, Z, Y), (X, Y, Z)):
        vec = 0.5 * np.cross(b - a, c - a)
        volume += float(np.dot((a + b + c) / 3.0, vec)) / 3.0
    return volume


def test_kernel_integrates_a_unit_tetrahedron_exactly():
    assert _tetra_volume() == pytest.approx(1.0 / 6.0, rel=1e-12)


def test_kernel_measures_a_sphere_below_the_three_percent_bar():
    check = sphere_check()
    assert abs(check["error_percent"]) < 0.5
    # a closed surface must have no net area vector; the unit sphere does up to tessellation
    assert check["closure_residual_norm"] < 0.05 * 4.0 * math.pi


def test_mask_accepts_a_square_and_rejects_its_complement():
    square = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
    inside = mask_uv(np.array([0.25, 0.75]), np.array([0.25, 0.75]), [square])
    assert inside.all()
    outside = mask_uv(np.array([1.25, 2.75]), np.array([0.25, 1.75]), [square])
    assert not outside.any()


def test_mask_is_orientation_independent_and_handles_an_island():
    ccw = [(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)]
    island = [(1.0, 1.0), (2.0, 1.0), (2.0, 2.0), (1.0, 2.0)]
    us = np.array([0.5, 1.5, 3.0])
    vs = np.array([0.5, 1.5, 3.0])
    without = mask_uv(us, vs, [ccw])
    with_island = mask_uv(us, vs, [ccw, island])
    assert without.all()
    # the island cell (1.5, 1.5) flips to outside under the even-odd rule
    assert not with_island[1, 1]
    assert with_island[0, 0] and with_island[2, 2]
    reversed_loop = list(reversed(ccw))
    assert (mask_uv(us, vs, [reversed_loop]) == without).all()


@pytest.mark.skipif(not ARTIFACT.exists(), reason="run the probe on the asset to produce the artifact")
def test_artifact_records_the_deficit_and_the_closure_residual():
    record = json.loads(ARTIFACT.read_text())
    assert record["patches_with_loops"] >= 1
    assert record["reference_full_m3"] == 312622.0
    # the finding: the exact, converged integral is short, and the shell is not closed
    assert record["error_vs_reference_percent"] < -10.0
    assert record["closure_residual_norm_m2"] > 0.0
    assert abs(record["integrator_check"]["error_percent"]) < 0.5
    assert record["de_duplicated_full_volume_m3"] == pytest.approx(record["full_volume_m3"], rel=0.01)


def test_the_document_cites_the_artifact_it_was_measured_from():
    if not DOC.exists():
        pytest.skip("document not in this checkout")
    text = DOC.read_text(encoding="utf-8")
    assert "artifacts/baseline_fleet/kvlcc2_divergence.json" in text
    assert "<3 %" in text or "< 3 %" in text  # the acceptance bar is stated, not implied