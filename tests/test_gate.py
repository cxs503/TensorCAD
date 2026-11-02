"""Tests for the hard-fail machine gate.

Every expected value here is external to the gate: the displaced volume comes
from :func:`tensorcad.hull.analytic_volume` (an integral of the analytic surface
that the mesh never saw), and the topology expectations (watertight, one
component, V-E+F=2) are geometric facts about a closed genus-0 shell, not
outputs of the code under test. The gate is exercised on deliberately broken
meshes as well, because a gate that has never rejected anything is not a gate.
"""
from __future__ import annotations

import math
from pathlib import Path

import pytest

from tensorcad import (
    HullParameters,
    Mesh,
    MeshGateError,
    analytic_volume,
    enforce_gate,
    gate_mesh,
    generate_hull,
)
from tensorcad import mesh as meshlib


def _hull(nx: int = 41, nz: int = 13) -> Mesh:
    return generate_hull(HullParameters(nx=nx, nz=nz)).mesh


def _exact(length: float = 10.0, beam: float = 2.0, draft: float = 0.8) -> float:
    return analytic_volume(length, beam, draft)


def _failed(report) -> set[str]:
    return {check.name for check in report.failures}


def test_the_generated_hull_passes_every_gate_check() -> None:
    """The reference model is the baseline: anything that fails here is a regression."""
    report = gate_mesh(_hull(), reference_volume=_exact())
    assert report.passed, report.to_dict()
    measurements = report.measurements
    assert measurements["boundary_edges"] == 0
    assert measurements["non_manifold_edges"] == 0
    assert measurements["components"] == 1
    assert measurements["euler_characteristic"] == 2
    assert measurements["zero_area_faces"] == 0
    # The analytic integral is the external reference, so the sign and size of the
    # measured deviation are meaningful: a coarse 41x13 grid under-integrates.
    assert -0.02 < measurements["relative_volume_error"] < 0.0


def test_gate_rejects_a_half_model_open_along_the_centreplane() -> None:
    """The historical failure: a starboard half exported as if it were a hull.

    It has faces, it renders, and it has no inside. Only a topology check finds it.
    """
    mesh = _hull()
    half = Mesh(
        vertices=mesh.vertices,
        faces=tuple(
            face
            for face in mesh.faces
            if meshlib.face_centroid(mesh, face)[1] >= 0.0
        ),
        units=mesh.units,
        coordinate_frame=mesh.coordinate_frame,
    )
    report = gate_mesh(half, reference_volume=_exact())
    assert not report.passed
    assert "watertight" in _failed(report)
    assert report.measurements["boundary_edges"] > 0
    with pytest.raises(MeshGateError, match="watertight"):
        enforce_gate(half, reference_volume=_exact(), label="half hull")


def test_gate_rejects_inverted_winding() -> None:
    mesh = _hull()
    inverted = Mesh(
        vertices=mesh.vertices,
        faces=tuple((a, c, b) for a, b, c in mesh.faces),
        units=mesh.units,
        coordinate_frame=mesh.coordinate_frame,
    )
    report = gate_mesh(inverted, reference_volume=_exact())
    assert not report.passed
    assert {"outward_orientation", "reference_volume"} <= _failed(report)
    assert report.measurements["signed_volume"] < 0.0


def test_gate_rejects_an_under_resolved_grid_against_the_analytic_volume() -> None:
    """A 178-triangle hull is closed and manifold -- and still 5.3 % wrong."""
    coarse = _hull(nx=11, nz=5)
    report = gate_mesh(coarse, reference_volume=_exact(), volume_tolerance=0.02)
    assert not report.passed
    assert "reference_volume" in _failed(report)
    detail = next(check.detail for check in report.failures if check.name == "reference_volume")
    assert "limit 2 %" in detail
    assert report.measurements["relative_volume_error"] < -0.05


def test_gate_enforces_an_explicit_facet_budget() -> None:
    """A resolvable mesh can still be too small for the consumer that asked."""
    report = gate_mesh(_hull(), min_faces=5000)
    assert not report.passed
    assert _failed(report) == {"facet_count"}


def test_gate_counts_components_and_the_euler_characteristic() -> None:
    mesh = _hull()
    shift = _exact()  # a second, disjoint hull translated clear of the first
    offset = 3.0 * max(abs(value) for point in mesh.vertices for value in point)
    two = Mesh(
        vertices=mesh.vertices
        + tuple((x + offset, y, z) for x, y, z in mesh.vertices),
        faces=mesh.faces
        + tuple((a + len(mesh.vertices), b + len(mesh.vertices), c + len(mesh.vertices))
                for a, b, c in mesh.faces),
        units=mesh.units,
        coordinate_frame=mesh.coordinate_frame,
    )
    report = gate_mesh(two, reference_volume=2.0 * shift)
    assert not report.passed
    assert {"component_count", "euler_characteristic"} <= _failed(report)
    assert report.measurements["components"] == 2
    assert report.measurements["euler_characteristic"] == 4
    accepted = gate_mesh(
        two, reference_volume=2.0 * shift, expected_components=2, expected_euler=4
    )
    assert accepted.passed, accepted.to_dict()


def test_gate_rejects_a_non_manifold_edge() -> None:
    mesh = _hull()
    extra = mesh.faces[0]
    duplicated = Mesh(
        vertices=mesh.vertices,
        faces=(*mesh.faces, extra),
        units=mesh.units,
        coordinate_frame=mesh.coordinate_frame,
    )
    report = gate_mesh(duplicated)
    assert not report.passed
    # "watertight" in validation terms means no boundary edge; the duplicated face
    # adds a *third* user of three edges, which the manifold check owns.
    assert "edge_manifold" in _failed(report)
    assert report.measurements["non_manifold_edges"] == 3
    assert not meshlib.is_watertight(duplicated)


def test_gate_rejects_a_zero_area_face() -> None:
    mesh = _hull()
    a, b, _ = mesh.faces[0]
    degenerate = Mesh(
        vertices=mesh.vertices,
        faces=(*mesh.faces, (a, b, a)),
        units=mesh.units,
        coordinate_frame=mesh.coordinate_frame,
    )
    report = gate_mesh(degenerate)
    assert not report.passed
    assert "no_degenerate_faces" in _failed(report)
    assert report.measurements["zero_area_faces"] >= 1


def test_gate_rejects_slivers_when_the_consumer_sets_a_limit() -> None:
    """The default hull ends in needle triangles; a caller that forbids them is told."""
    report = gate_mesh(_hull(), max_aspect=4.0)
    assert not report.passed
    assert _failed(report) == {"facet_aspect"}
    assert report.measurements["max_facet_aspect"] > 4.0


def test_gate_rejects_non_finite_reference_values() -> None:
    for bad in (0.0, -1.0, math.nan, math.inf):
        with pytest.raises(ValueError, match="reference_volume"):
            gate_mesh(_hull(), reference_volume=bad)


def test_gate_rejects_impossible_expectations() -> None:
    with pytest.raises(ValueError, match="expected_components"):
        gate_mesh(_hull(), expected_components=0)
    with pytest.raises(ValueError, match="volume_tolerance"):
        gate_mesh(_hull(), volume_tolerance=-1.0)
    with pytest.raises(ValueError, match="min_faces"):
        gate_mesh(_hull(), min_faces=0)


def test_report_serialises_without_losing_a_check(tmp_path: Path) -> None:
    import json

    report = gate_mesh(_hull(), reference_volume=_exact())
    payload = json.loads(json.dumps(report.to_dict()))
    assert payload["schema"] == "tensorcad.mesh.gate/1"
    assert payload["passed"] is True
    assert len(payload["checks"]) == len(report.checks)
    assert payload["measurements"]["euler_characteristic"] == 2
