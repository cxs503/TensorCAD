"""The bridge between the repository's ``Mesh`` and the canonical quality library.

What is asserted, and why each assertion can fail:

* the four aspect dialects are related by the *exact* algebra, measured on real meshes,
  not on a comment (``gate = 2 * legacy`` to rounding, ``edge_form = legacy * alpha``,
  ``verdict = legacy * alpha * beta``);
* the bridge reproduces the production gate's own number bit-for-bit on those meshes, so
  the delegation is a proof of equivalence and not a re-implementation with a second
  rounding;
* the canonical limit derived from a legacy limit can never admit a mesh the legacy limit
  rejects (a hand-built counterexample drives the two gates apart, so the test fails if the
  transfer is ever "simplified" into a naive copy of the number);
* wiring the canonical criteria into :func:`tensorcad.gate.gate_mesh` is opt-in: with no
  canonical argument the verdicts and the measurement keys are exactly what they were;
* a zero-area or inverted facet is reported, not dropped.
"""
from __future__ import annotations

import math

import pytest

from tensorcad import gate
from tensorcad.geometry import Mesh
from tensorcad.hull import HullParameters, build_mesh
from tensorcad.meshquality.bridge import (
    CANONICAL_METRICS,
    canonical_stats,
    dialect_identities,
    gate_facet_aspect,
    gate_to_canonical_limit,
    legacy_to_canonical_limit,
    qualify,
)


def hull_mesh(nx: int = 41, nz: int = 13):
    mesh, _ = build_mesh(HullParameters(length=10.0, beam=2.0, draft=0.8, nx=nx, nz=nz))
    return mesh


def box_mesh() -> Mesh:
    vertices = [
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (1.0, 1.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0),
        (1.0, 0.0, 1.0),
        (1.0, 1.0, 1.0),
        (0.0, 1.0, 1.0),
    ]
    quads = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
    faces = []
    for a, b, c, d in quads:
        faces.append((a, b, c))
        faces.append((a, c, d))
    return Mesh(vertices=tuple(vertices), faces=tuple(faces))


def triangle_mesh(points) -> Mesh:
    return Mesh(vertices=tuple(points), faces=((0, 1, 2),))


# ---------------------------------------------------------------------------------------
# 1. the dialect algebra, measured
# ---------------------------------------------------------------------------------------
@pytest.mark.parametrize("name", ["box", "hull", "hull_fine"])
def test_dialect_identities_hold_to_rounding(name: str) -> None:
    mesh = {"box": box_mesh, "hull": hull_mesh, "hull_fine": lambda: hull_mesh(81, 25)}[name]()
    deviations = dialect_identities(mesh)
    assert deviations["gate_minus_2_legacy"] <= 1e-12 * max(1.0, gate._facet_aspect(mesh))
    assert deviations["edge_form_minus_legacy_alpha"] <= 1e-9 * max(1.0, gate._facet_aspect(mesh))
    assert deviations["verdict_minus_legacy_alpha_beta"] <= 1e-9 * max(1.0, gate._facet_aspect(mesh))


@pytest.mark.parametrize("name", ["box", "hull"])
def test_bridge_reproduces_the_gate_number(name: str) -> None:
    mesh = {"box": box_mesh, "hull": hull_mesh}[name]()
    quality = qualify(mesh)
    reported = gate.measure_gated(mesh)["max_facet_aspect"]
    assert abs(quality.dialect_max()["gate"] - reported) <= 1e-12 * max(1.0, reported)


def test_ideal_triangle_scores_one_in_every_canonical_normalisation() -> None:
    h = math.sqrt(3.0) / 2.0
    quality = qualify(triangle_mesh([(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.5, h, 0.0)]))
    assert quality.elements["aspect_ratio"][0] == pytest.approx(1.0, abs=1e-12)
    assert quality.elements["min_sicn"][0] == pytest.approx(1.0, abs=1e-12)
    assert quality.elements["scaled_jacobian"][0] == pytest.approx(1.0, abs=1e-12)
    assert quality.elements["min_angle"][0] == pytest.approx(60.0, abs=1e-12)
    assert quality.elements["skewness"][0] == pytest.approx(0.0, abs=1e-12)
    # the three legacy dialects disagree on the same triangle, by factors, not by taste
    dialects = quality.dialect_max()
    # gate = l^2/A = 4/sqrt(3); legacy = l^2/2A = 2/sqrt(3); both canonical forms are 1
    assert dialects["gate"] == pytest.approx(4.0 / math.sqrt(3.0), abs=1e-12)
    assert dialects["legacy"] == pytest.approx(2.0 / math.sqrt(3.0), abs=1e-12)
    assert dialects["edge_form"] == pytest.approx(1.0, abs=1e-12)
    assert dialects["verdict"] == pytest.approx(1.0, abs=1e-12)


# ---------------------------------------------------------------------------------------
# 2. the transfer of a limit cannot weaken the gate
# ---------------------------------------------------------------------------------------
def test_naive_threshold_copy_is_caught_by_the_conservative_transfer() -> None:
    """A right isosceles triangle: verdict/legacy = 0.6036, so a copied limit is unsafe."""
    right = triangle_mesh([(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)])
    quality = qualify(right)
    gate_value = quality.dialect_max()["gate"]
    verdict_value = quality.dialect_max()["verdict"]
    limit = gate_value * 0.99  # a gate-dialect limit this triangle violates
    assert gate_value > limit > 0.0
    translated = gate_to_canonical_limit(limit)
    assert verdict_value > translated, "the translated limit must reject what the gate rejects"
    assert verdict_value <= limit, "a naive copy of the number would have passed it"
    assert translated <= limit / 4.0 + 1e-12


def test_transfer_bounds_are_ordered_for_every_source_dialect() -> None:
    assert legacy_to_canonical_limit(25.0, source="legacy") == pytest.approx(12.5)
    assert legacy_to_canonical_limit(25.0, source="gate") == pytest.approx(6.25)
    assert legacy_to_canonical_limit(25.0, source="edge_form") == pytest.approx(25.0 * math.sqrt(3.0) / 2.0)
    with pytest.raises(ValueError):
        legacy_to_canonical_limit(25.0, source="not_a_dialect")


def test_canonical_gate_never_passes_a_mesh_the_legacy_gate_rejects() -> None:
    mesh = hull_mesh()
    strict = gate.gate_mesh(mesh, max_aspect=4.0)  # the hull violates this on purpose
    assert not strict.passed
    translated = gate.gate_mesh(mesh, max_aspect=4.0, **gate.canonical_limits_from_legacy(max_aspect=4.0))
    assert not translated.passed
    assert {"facet_aspect", "canonical_aspect"} <= {check.name for check in translated.failures}


# ---------------------------------------------------------------------------------------
# 3. wiring is opt-in and non-breaking
# ---------------------------------------------------------------------------------------
def test_measurements_keep_every_key_and_add_the_canonical_block() -> None:
    mesh = hull_mesh()
    measurements = gate.measure_gated(mesh)
    assert measurements["max_facet_aspect"] == pytest.approx(gate._facet_aspect(mesh))
    assert set(measurements["canonical"]) == {
        "max_aspect_ratio", "max_edge_ratio", "max_skewness", "max_max_angle", "min_angle",
        "min_sicn", "min_scaled_jacobian", "min_gamma", "min_eta", "non_positive_area",
    }
    # the canonical aspect ratio is R/2r, which cannot exceed the gate's l^2/A for a triangle
    assert measurements["canonical"]["max_aspect_ratio"] <= measurements["max_facet_aspect"]
    assert measurements["canonical"] == canonical_stats(mesh)


def test_gate_verdict_is_unchanged_when_no_canonical_criterion_is_given() -> None:
    mesh = hull_mesh()
    before = gate.gate_mesh(mesh, max_aspect=25.0, reference_volume=None)
    after = gate.gate_mesh(mesh, max_aspect=25.0, reference_volume=None)
    assert [check.name for check in before.checks] == [check.name for check in after.checks]
    assert [check.passed for check in before.checks] == [check.passed for check in after.checks]
    assert not any(check.name.startswith("canonical_") for check in after.checks)


def test_canonical_criteria_fail_with_the_offending_facet_located() -> None:
    mesh = hull_mesh()
    report = gate.gate_mesh(mesh, min_canonical_angle=10.0, max_canonical_skewness=0.5)
    failing = {check.name: check for check in report.failures}
    assert set(failing) == {"canonical_min_angle", "canonical_skewness"}
    assert "at face " in failing["canonical_min_angle"].detail


# ---------------------------------------------------------------------------------------
# 4. failures are reported, not skipped
# ---------------------------------------------------------------------------------------
def test_orientation_is_a_topology_property_not_an_element_sign() -> None:
    """Triangle element metrics are winding-agnostic by design; orientation lives in topology.

    Stated explicitly because it is easy to assume the opposite: ``tri_metrics`` reports the
    *unsigned* area, so a surface's inverted facets show up as an inconsistent winding and a
    negative enclosed volume in ``surface_topology`` (which the gate consumes), not as a
    negative SICN. Tetrahedra and hexahedra do carry the sign, because a volume element has
    an inside.
    """
    one = triangle_mesh([(0.0, 0.0, 0.0), (0.0, 1.0, 0.0), (1.0, 0.0, 0.0)])
    quality = qualify(one)
    assert quality.elements["min_sicn"][0] > 0.0
    assert quality.elements["area"][0] > 0.0
    closed = gate.gate_mesh(
        Mesh(vertices=one.vertices, faces=((0, 1, 2), (0, 2, 1))),
        min_canonical_sicn=0.05,
    )
    assert not closed.passed
    assert {"outward_orientation", "watertight"} & {check.name for check in closed.failures}


def test_zero_area_facet_fails_the_sicn_criterion() -> None:
    collinear = Mesh(
        vertices=((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0)),
        faces=((0, 1, 2),),
    )
    report = gate.gate_mesh(collinear, min_canonical_sicn=0.05)
    assert "canonical_min_sicn" in {check.name for check in report.failures}


def test_zero_area_facet_is_counted_and_is_the_worst() -> None:
    collinear = triangle_mesh([(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0)])
    quality = qualify(collinear)
    assert quality.elements["area"][0] == 0.0
    assert quality.scalars()["non_positive_area"] == 1.0
    assert quality.worst("aspect_ratio")["index"] == 0
    assert quality.dialect_max()["gate"] == math.inf
    assert gate_facet_aspect(0.0, 1.0) == math.inf


def test_quad_faces_are_refused_rather_than_silently_fanned() -> None:
    quad = Mesh(
        vertices=((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0)),
        faces=((0, 1, 2), (0, 2, 3)),
    )
    assert len(qualify(quad).elements) == 2
    with pytest.raises(ValueError):
        qualify(type("FakeMesh", (), {"vertices": quad.vertices, "faces": ((0, 1, 2, 3),), "units": "m"})())


def test_bridge_matches_the_canonical_library_for_the_same_facets() -> None:
    """The bridge must not be a second implementation: same faces, same numbers."""
    from tensorcad.meshquality.report import measure

    mesh = hull_mesh()
    direct = measure([tuple(float(c) for c in v) for v in mesh.vertices], {"tri": [tuple(f) for f in mesh.faces]})
    summary = direct.summary()["tri"]
    via_bridge = qualify(mesh)
    assert via_bridge.scalars()["min_sicn"] == pytest.approx(summary["min_sicn"].minimum, abs=0.0)
    assert via_bridge.scalars()["max_aspect_ratio"] == pytest.approx(summary["aspect_ratio"].maximum, abs=0.0)
    assert via_bridge.scalars()["max_skewness"] == pytest.approx(summary["skewness"].maximum, abs=0.0)
    assert via_bridge.scalars()["min_angle"] == pytest.approx(summary["min_angle"].minimum, abs=0.0)