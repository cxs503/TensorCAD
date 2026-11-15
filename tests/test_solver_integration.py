"""Tests for the unified solver-interface integration (:mod:`tensorcad.solver_integration`).

What these pin:
- all three routes (FVM, DG, LBM) declare the same z-up frame on the same geometry;
- the analytic volume of the ellipsoid is the reference every route is scored against;
- a y-up tamper is refused by every route's verifier (the 5415 lesson);
- the integration report is self-consistent and machine-readable.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from tensorcad.solver_integration import (
    INTEGRATION_SCHEMA,
    IntegrationReport,
    RouteResult,
    frame_tamper_refused,
    run_integration,
    summary_markdown,
    verify_frame_consistency,
    _analytic_volume,
)
from tensorcad.solver_interface import COORDINATE_CONTRACT, HandoffError


# --------------------------------------------------------------------------- analytic reference


def test_analytic_volume_is_4_3_pi_abc():
    """4/3 π a b c — the reference every route is scored against."""
    assert _analytic_volume((1.0, 1.0, 1.0)) == pytest.approx(4.0 / 3.0 * math.pi)
    assert _analytic_volume((2.0, 1.0, 0.5)) == pytest.approx(4.0 / 3.0 * math.pi * 2.0 * 1.0 * 0.5)


# --------------------------------------------------------------------------- frame consistency


def test_frame_consistency_passes_when_all_zup():
    routes = [
        RouteResult(route="volume_msh", status="ran", metrics={"volume_rel_error": 0.01},
                    frame="z-up"),
        RouteResult(route="curved_wall", status="ran", metrics={"circumference_rel_error": 1e-14},
                    frame="z-up"),
        RouteResult(route="lattice", status="ran", metrics={"area_rel_error": 0.05},
                    frame="z-up"),
    ]
    assert verify_frame_consistency(routes) is True


def test_frame_consistency_fails_on_yup():
    routes = [
        RouteResult(route="volume_msh", status="ran", metrics={"volume_rel_error": 0.01},
                    frame="z-up"),
        RouteResult(route="curved_wall", status="ran", metrics={"circumference_rel_error": 1e-14},
                    frame="y-up"),
    ]
    assert verify_frame_consistency(routes) is False


def test_frame_consistency_fails_on_bad_volume():
    routes = [
        RouteResult(route="volume_msh", status="ran", metrics={"volume_rel_error": 0.10},
                    frame="z-up"),
    ]
    assert verify_frame_consistency(routes) is False


def test_frame_consistency_fails_on_bad_circumference():
    routes = [
        RouteResult(route="curved_wall", status="ran",
                    metrics={"circumference_rel_error": 0.01}, frame="z-up"),
    ]
    assert verify_frame_consistency(routes) is False


# --------------------------------------------------------------------------- frame tamper


def test_frame_tamper_refused_volume(tmp_path: Path):
    """A y-up manifest on the volume route is refused."""
    route_dir = tmp_path / "volume_msh"
    route_dir.mkdir()
    assert frame_tamper_refused(route_dir, "volume_msh") is True


def test_frame_tamper_refused_lattice(tmp_path: Path):
    """A y-up manifest on the lattice route is refused."""
    route_dir = tmp_path / "lattice"
    route_dir.mkdir()
    assert frame_tamper_refused(route_dir, "lattice") is True


def test_frame_tamper_refused_dg(tmp_path: Path):
    """A y-up manifest on the DG route is refused."""
    route_dir = tmp_path / "dg"
    route_dir.mkdir()
    assert frame_tamper_refused(route_dir, "curved_wall") is True


# --------------------------------------------------------------------------- full integration


@pytest.mark.parametrize("radial_cells", [6, 10])
def test_integration_all_routes_pass(tmp_path: Path, radial_cells: int):
    """All three routes run on the same ellipsoid and pass frame checks."""
    out = tmp_path / "integration"
    report = run_integration(out, semi_axes=(1.0, 0.5, 0.3),
                             fvm_resolution=(18, 24),
                             lbm_radial_cells=radial_cells)

    assert report.schema == INTEGRATION_SCHEMA
    assert report.semi_axes == (1.0, 0.5, 0.3)
    assert report.analytic_volume_m3 == pytest.approx(_analytic_volume((1.0, 0.5, 0.3)))
    assert len(report.routes) == 3
    assert report.frame_consistent is True
    assert report.frame_tamper_refused is True
    assert report.all_passed is True

    # Check each route
    route_map = {r.route: r for r in report.routes}
    assert route_map["volume_msh"].status == "ran"
    assert route_map["volume_msh"].frame == "z-up"
    assert route_map["volume_msh"].metrics["volume_rel_error"] < 0.05

    assert route_map["curved_wall"].status == "ran"
    assert route_map["curved_wall"].frame == "z-up"
    assert route_map["curved_wall"].metrics["circumference_rel_error"] < 1e-10

    assert route_map["lattice"].status == "ran"
    assert route_map["lattice"].frame == "z-up"
    assert route_map["lattice"].metrics["area_rel_error"] < 0.10


def test_integration_writes_artifacts(tmp_path: Path):
    """The integration run writes a JSON report and a markdown summary."""
    out = tmp_path / "integration"
    run_integration(out, semi_axes=(0.5, 0.25, 0.15),
                    fvm_resolution=(12, 16), lbm_radial_cells=6)

    report_json = out / "integration_report.json"
    summary_md = out / "summary.md"
    assert report_json.exists()
    assert summary_md.exists()

    payload = json.loads(report_json.read_text())
    assert payload["schema"] == INTEGRATION_SCHEMA
    assert payload["all_passed"] is True
    assert payload["frame_consistent"] is True
    assert payload["frame_tamper_refused"] is True

    md = summary_md.read_text()
    assert "PASS" in md
    assert "z-up" in md
    assert "5415" in md


def test_summary_markdown_contains_key_sections():
    report = IntegrationReport(
        schema=INTEGRATION_SCHEMA,
        semi_axes=(1.0, 0.5, 0.3),
        analytic_volume_m3=0.628319,
        routes=(
            RouteResult(route="volume_msh", status="ran",
                        metrics={"signed_volume_m3": 0.628, "volume_rel_error": 0.0005},
                        frame="z-up"),
            RouteResult(route="curved_wall", status="ran",
                        metrics={"circumference_m": 6.283, "circumference_rel_error": 1e-14},
                        frame="z-up"),
            RouteResult(route="lattice", status="ran",
                        metrics={"cross_section_area_m2": 3.14, "area_rel_error": 0.001},
                        frame="z-up"),
        ),
        frame_consistent=True,
        frame_tamper_refused=True,
        all_passed=True,
        elapsed_s=1.5,
    )
    md = summary_markdown(report)
    assert "Solver Interface Integration Report" in md
    assert "volume_msh" in md
    assert "curved_wall" in md
    assert "lattice" in md
    assert "Frame Consistency" in md
    assert "5415 Lesson" in md
    assert "PASS" in md


def test_integration_report_is_serialisable(tmp_path: Path):
    out = tmp_path / "integration"
    report = run_integration(out, semi_axes=(0.3, 0.2, 0.1),
                             fvm_resolution=(12, 16), lbm_radial_cells=6)
    # Round-trip through JSON
    payload = report.to_dict()
    text = json.dumps(payload)
    loaded = json.loads(text)
    assert loaded["schema"] == INTEGRATION_SCHEMA
    assert loaded["all_passed"] is True
