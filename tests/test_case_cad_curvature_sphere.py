"""CASE-CAD-CURVATURE-SPHERE: the discrete mean curvature of a sphere against ``2/R``.

The verified case owns the full series (``benchmarks/verified/sphere_mean_curvature``);
this test pins the two claims that must not silently drift -- the integral reading is
inside the 2 % criterion at a practical resolution, and the pointwise reading is second
order where the triangulation is structured -- at a size that runs in seconds.
"""
from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parents[1]
_CASE = _REPO_ROOT / "benchmarks" / "verified" / "sphere_mean_curvature" / "run_benchmark.py"


def _load_case():
    spec = importlib.util.spec_from_file_location("case_sphere_mean_curvature", _CASE)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def case():
    return _load_case()


@pytest.fixture(scope="module")
def result(case):
    """One quick run, shared by the assertions below (~15 s)."""
    return case.run(
        radius=1.0,
        levels=(16, 32),
        condition_levels=(32,),
        resolutions=((8, 16), (16, 32), (32, 64)),
        radii=(0.5, 1.0),
        radius_level=16,
        probes=24,
    )


def test_the_analytic_reference_is_the_sum_of_the_principal_curvatures(case) -> None:
    """2/R, and its surface integral is 8 pi R -- the case's own error bar."""
    assert case.reference_sum(2.0) == pytest.approx(1.0, rel=0.0)
    assert case.reference_total(1.0) == pytest.approx(8.0 * math.pi, rel=1e-15)
    assert case.reference_total(3.0) == pytest.approx(24.0 * math.pi, rel=1e-15)


def test_the_integral_reading_is_inside_the_two_percent_criterion(case, result) -> None:
    headline = result["headline"]
    assert headline["reference"] == pytest.approx(2.0, rel=0.0)
    assert abs(headline["relative_error"]) < 0.02, headline
    # and it is not merely inside: on a closed mesh the cotangent divergence
    # telescopes, so the integral is a fraction of a percent even at level 16
    assert abs(headline["relative_error"]) < 5e-3, headline


def test_the_pointwise_reading_needs_a_well_shaped_triangulation(case, result) -> None:
    """Structured: inside 2 %. Raw marching tets: an order of magnitude outside."""
    structured = result["structured"]["series"]
    fine = next(entry for entry in structured if (entry["nu"], entry["nv"]) == (16, 32))
    assert fine["pointwise_median_relative_error"] < 0.02, fine
    assert result["headline"]["pointwise_median_relative_error_raw"] > 0.2
    assert result["headline"]["pointwise_median_relative_error_conditioned"] < 0.03


def test_the_pointwise_estimate_keeps_second_order_on_a_structured_sphere(result) -> None:
    order = result["structured"]["observed_order"]
    assert order == pytest.approx(2.0, abs=0.5), order


def test_the_relative_error_does_not_depend_on_the_radius(case, result) -> None:
    series = result["radius"]["series"]
    errors = [entry["integral_relative_error"] for entry in series]
    assert max(errors) - min(errors) < 1e-9, errors
    for entry in series:
        assert abs(entry["integral_relative_error"]) < 0.02, entry


def test_the_field_path_is_the_reference_for_the_measurement_chain(result) -> None:
    """Every G4 row: the analytic estimator returns 2/R to round-off."""
    field_rows = [row for row in result["reference_chain"] if row["gate"] == "G4"]
    assert len(field_rows) == 7
    for row in field_rows:
        assert abs(row["error"]) < 1e-6, row


def test_every_gate_of_the_case_passes(result) -> None:
    assert result["failing"] == []
    assert result["all_gates_pass"] is True
    # the raw-pointwise diagnostic is expected to break the criterion; it is reported
    assert len(result["diagnostics"]) == 1
    assert result["diagnostics"][0]["holds"] is False


def test_the_uv_sphere_fixture_is_a_legal_closed_mesh(case) -> None:
    """No zero-area face, no collapsed pole: otherwise the fixture measures itself."""
    from tensorcad import mesh as meshlib

    mesh = case.uv_sphere(1.0, 16, 32)
    assert len(mesh.vertices) == (16 - 1) * 32 + 2  # (nu-1) rings of nv plus two poles
    assert list(meshlib.zero_area_faces(mesh)) == []
    result = case.measure(mesh, 1.0)
    assert result["area_relative_error"] == pytest.approx(0.0, abs=1e-2)
    assert result["pointwise"]["count"] == len(mesh.vertices)