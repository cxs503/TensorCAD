"""Gold-regression and reference-integrity tests for GAP-CAD-T1-IGES-PIPE.

The first published version of the benchmark judged the IGES import pipeline against a
5x5 Gauss quadrature on the *whole* patch domain.  That reference under-integrates the
keel singularity of the analytic baselines (``y ~ zeta**(1/m)``, ``m`` up to 24) by
3.5 %, i.e. the reference was the defect: Series 60's exact block coefficient gives
V = 3.75 m3 for Cb = 0.60, the old reference said 3.61999 m3, and the pipeline was then
"measured" as converging to the wrong number.

These tests pin the repaired reference:
* it reproduces the analytic closed forms of the baselines (the independent anchor);
* its own panel ladder settles, so a converged value is a number, not a rung;
* the pipeline, at the density the benchmark uses, stays inside the 3 % gate of the
  reference of the file it actually read.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "benchmarks" / "verified" / "iges_import"))

from gold_loft import (  # noqa: E402
    PatchDerivatives,
    brep_quadrature,
    graded_panels,
    loft_gold_hull,
    read_stl,
)

from tensorcad.geomquality.baselines import kvlcc2_class, series60  # noqa: E402
from tensorcad.geometry import load_cad_geometry  # noqa: E402
from tensorcad.geomquality.io import write_iges  # noqa: E402

GOLD = Path("/root/ship-performance-platform/backend/data/geometry/ships/"
            "DTMB5415_Hull_exact.stl")


def _volume(vertices, faces) -> float:
    points = np.asarray(vertices, dtype=float)
    triangles = np.asarray(faces, dtype=int)
    a, b, c = points[triangles[:, 0]], points[triangles[:, 1]], points[triangles[:, 2]]
    return float(np.einsum("ij,ij->i", a, np.cross(b, c)).sum() / 6.0)


# --------------------------------------------------------------------------- #
# the cached derivative evaluator must be a cache, not a reimplementation
# --------------------------------------------------------------------------- #
def test_cached_derivatives_equal_the_library():
    surface = series60(0.60).brep.faces[0].surface
    cached = PatchDerivatives(surface)
    for u, v in ((0.0, 0.0), (0.13, 0.77), (0.5, 0.5), (0.91, 0.02), (1.0, 1.0)):
        point, du, dv = cached(u, v)
        reference = surface.derivatives(u, v, order=1)
        assert np.allclose(point, reference[(0, 0)], atol=1e-12)
        assert np.allclose(du, reference[(1, 0)], atol=1e-12)
        assert np.allclose(dv, reference[(0, 1)], atol=1e-12)


# --------------------------------------------------------------------------- #
# graded panel grid
# --------------------------------------------------------------------------- #
def test_graded_panels_tile_the_domain_and_cluster_at_both_ends():
    for count in (1, 2, 3, 4, 8, 9):
        panels = graded_panels(-2.0, 3.0, count, ratio=0.5)
        assert len(panels) == count
        assert panels[0][0] == pytest.approx(-2.0)
        assert panels[-1][1] == pytest.approx(3.0)
        for (a, b), (c, d) in zip(panels, panels[1:]):
            assert b == pytest.approx(c)
            assert b > a
        if count >= 4:
            widths = [b - a for a, b in panels]
            assert widths[0] < widths[count // 2]   # clustered at the low end
            assert widths[-1] < widths[count // 2]  # clustered at the high end


# --------------------------------------------------------------------------- #
# the reference integral reproduces the analytic volumes
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("specimen", [series60(0.60), series60(0.65), kvlcc2_class()],
                         ids=["series60_cb060", "series60_cb065", "kvlcc2"])
def test_reference_matches_the_analytic_closed_form(specimen):
    closed = specimen.reference["closed_form"]["volume"]
    volume = brep_quadrature(specimen.brep, order=5, panels=4)["volume"]
    assert abs(volume - closed) / closed < 5e-3, (volume, closed)


def test_reference_ladder_settles_and_the_old_single_panel_does_not():
    specimen = series60(0.60)
    closed = specimen.reference["closed_form"]["volume"]
    rungs = [brep_quadrature(specimen.brep, panels=p)["volume"] for p in (1, 2, 4)]
    # the previous reference (one panel over the whole domain) is off by >3 %: that is
    # exactly the defect this test exists for.
    assert abs(rungs[0] - closed) / closed > 0.03
    # the graded ladder converges to the analytic value.
    assert abs(rungs[-1] - closed) / closed < 5e-3
    assert abs(rungs[-1] - rungs[-2]) < abs(rungs[1] - rungs[0])


# --------------------------------------------------------------------------- #
# end-to-end: pipeline vs the reference of the file it read
# --------------------------------------------------------------------------- #
def test_pipeline_reproduces_the_written_patch_volume_within_3pct(tmp_path):
    specimen = series60(0.60)
    target = tmp_path / "series60_cb060.iges"
    write_iges(target, specimen.brep)
    reference = brep_quadrature(specimen.brep, order=5, panels=4)["volume"]
    mesh, report = load_cad_geometry(target, u_samples=96, v_samples=48)
    assert report["after_repair"]["watertight"]
    # the report publishes ``components``; the benchmark driver once read
    # ``component_count`` and died with KeyError before recording a single row.
    assert report["after_repair"]["components"] == 1
    assert report["after_repair"]["closed_components"] == 1
    measured = _volume(mesh.vertices, mesh.faces)
    assert abs(measured - reference) / reference < 0.03, (measured, reference)


@pytest.mark.skipif(not GOLD.exists(), reason="platform gold mesh not present")
def test_dtmb5415_loft_imports_watertight_and_tracks_the_gold_surface(tmp_path):
    triangles = read_stl(GOLD)
    brep, info = loft_gold_hull(triangles)
    target = tmp_path / "gold_loft.iges"
    write_iges(target, brep)
    mesh, report = load_cad_geometry(target, u_samples=64, v_samples=32)
    assert report["after_repair"]["watertight"]
    assert report["after_repair"]["components"] == 1
    # The loft is an interpolation of the gold sections, so the imported body must sit on
    # the gold surface: bounding box agreement to 1 % of the length.
    points = np.asarray(mesh.vertices)
    length = float(np.ptp(triangles[:, :, 0]))
    assert abs(np.ptp(points[:, 0]) - length) / length < 0.01
    assert info["same_sense"] and len(info["same_sense"]) == 4


def test_benchmark_driver_no_longer_reads_the_absent_key():
    source = (ROOT / "benchmarks" / "verified" / "iges_import" / "run_benchmark.py").read_text()
    assert '"component_count"' not in source
    assert '"components"' in source
