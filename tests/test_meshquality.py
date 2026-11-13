"""Tests for the mesh-quality library: definitions, invariance, gate and references.

The assertions are mostly *exact* values on canonical elements, because that is what
distinguishes a metric library from a number printer: an equilateral triangle has
``min_sicn == 1`` by the definition in ``tensorcad.meshquality.metrics``, and a regular
tetrahedron has ``gmsh_icn == (1/sqrt(2))**(2/3)``. A drift in a normalisation shows up
here as a fraction of a percent, which a "looks reasonable" test would miss.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from tensorcad.meshquality import cli, reference
from tensorcad.meshquality.io import MeshSpec, load_spec, write_spec
from tensorcad.meshquality.metrics import (
    _corner_quality,
    _hex_jacobian,
    hex_metrics,
    quad_metrics,
    tet_metrics,
    tri_metrics,
)
from tensorcad.meshquality.report import (
    DEFAULT_GATE,
    MetricLimit,
    QualityGate,
    measure,
    summary_stats,
)


# --------------------------------------------------------------------------------------
# canonical elements
# --------------------------------------------------------------------------------------
def regular_tet_points() -> list[tuple[float, float, float]]:
    h = math.sqrt(3.0) / 2.0
    return [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.5, h, 0.0), (0.5, h / 3.0, math.sqrt(2.0 / 3.0))]


def test_equilateral_triangle_scores_one() -> None:
    points = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.5, math.sqrt(3.0) / 2.0, 0.0)]
    quality = tri_metrics(points, [(0, 1, 2)])
    assert quality["min_sicn"][0] == pytest.approx(1.0, abs=1e-12)
    assert quality["scaled_jacobian"][0] == pytest.approx(1.0, abs=1e-12)
    assert quality["gamma"][0] == pytest.approx(1.0, abs=1e-12)
    assert quality["aspect_ratio"][0] == pytest.approx(1.0, abs=1e-12)
    assert quality["eta"][0] == pytest.approx(1.0, abs=1e-12)
    assert quality["min_angle"][0] == pytest.approx(60.0, abs=1e-12)
    assert quality["area"][0] == pytest.approx(math.sqrt(3.0) / 4.0, rel=1e-15)


def test_right_isosceles_triangle_scores_known_values() -> None:
    # Half of a unit square: area 1/2, angles 45/45/90, 2r/R = 2 (sqrt(2) - 1).
    points = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)]
    quality = tri_metrics(points, [(0, 1, 2)])
    assert quality["min_angle"][0] == pytest.approx(45.0, abs=1e-12)
    assert quality["area"][0] == pytest.approx(0.5, rel=1e-15)
    assert quality["gamma"][0] == pytest.approx(2.0 * (math.sqrt(2.0) - 1.0), abs=1e-12)
    assert quality["scaled_jacobian"][0] == pytest.approx(math.sqrt(2.0) / math.sqrt(3.0), abs=1e-12)
    # Against the equilateral reference the right triangle is a real distortion:
    # 2 det(M)/||M||^2 with det(M) = 4A/sqrt(3) = 2/sqrt(3), ||M||^2 = 8/3.
    assert quality["min_sicn"][0] == pytest.approx(math.sqrt(3.0) / 2.0, abs=1e-12)


def test_three_four_five_triangle() -> None:
    points = [(0.0, 0.0, 0.0), (4.0, 0.0, 0.0), (0.0, 3.0, 0.0)]
    quality = tri_metrics(points, [(0, 1, 2)])
    assert quality["area"][0] == pytest.approx(6.0, rel=1e-15)
    assert quality["edge_ratio"][0] == pytest.approx(5.0 / 3.0, rel=1e-15)
    # R = 2.5 (hypotenuse/2), r = 1 (area/s), so Verdict aspect ratio = 2.5 / 2 = 1.25.
    assert quality["aspect_ratio"][0] == pytest.approx(1.25, rel=1e-14)
    assert quality["gamma"][0] == pytest.approx(0.8, rel=1e-14)
    assert quality["min_angle"][0] == pytest.approx(math.degrees(math.atan2(3.0, 4.0)), abs=1e-12)


def test_unit_square_quad_scores_one() -> None:
    points = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0)]
    quality = quad_metrics(points, [(0, 1, 2, 3)])
    assert quality["min_sicn"][0] == pytest.approx(1.0, abs=1e-12)
    assert quality["scaled_jacobian"][0] == pytest.approx(1.0, abs=1e-12)
    assert quality["eta"][0] == pytest.approx(1.0, abs=1e-12)
    assert quality["min_angle"][0] == pytest.approx(90.0, abs=1e-12)
    assert quality["area"][0] == pytest.approx(1.0, rel=1e-15)


def test_stretched_quad_keeps_scaled_jacobian_but_loses_sicn() -> None:
    # 8x1 rectangle: orthogonal corners (shape score 1) but the condition number 2*8*1/(64+1)=16/65.
    points = [(0.0, 0.0, 0.0), (8.0, 0.0, 0.0), (8.0, 1.0, 0.0), (0.0, 1.0, 0.0)]
    quality = quad_metrics(points, [(0, 1, 2, 3)])
    assert quality["scaled_jacobian"][0] == pytest.approx(1.0, abs=1e-12)
    assert quality["min_sicn"][0] == pytest.approx(16.0 / 65.0, rel=1e-12)


def test_unit_cube_hex_scores_one() -> None:
    points = [
        (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 1.0), (0.0, 1.0, 1.0),
    ]
    quality = hex_metrics(points, [tuple(range(8))])
    assert quality["volume"][0] == pytest.approx(1.0, rel=1e-15)
    assert quality["min_sicn"][0] == pytest.approx(1.0, abs=1e-12)
    assert quality["scaled_jacobian"][0] == pytest.approx(1.0, abs=1e-12)
    assert quality["min_angle"][0] == pytest.approx(90.0, abs=1e-12)
    assert quality["max_angle"][0] == pytest.approx(90.0, abs=1e-12)
    assert quality["edge_ratio"][0] == pytest.approx(1.0, rel=1e-15)


def test_regular_tetrahedron_scores_one_in_the_sicn_unit() -> None:
    quality = tet_metrics(regular_tet_points(), [(0, 1, 2, 3)])
    assert quality["volume"][0] == pytest.approx(math.sqrt(2.0) / 12.0, rel=1e-14)
    assert quality["min_sicn"][0] == pytest.approx(1.0, abs=1e-12)
    assert quality["gamma"][0] == pytest.approx(1.0, abs=1e-12)
    assert quality["eta"][0] == pytest.approx(1.0, abs=1e-12)
    assert quality["scaled_jacobian"][0] == pytest.approx(1.0 / math.sqrt(2.0), rel=1e-12)
    assert quality["min_dihedral"][0] == pytest.approx(math.degrees(math.acos(1.0 / 3.0)), abs=1e-12)


def test_gmsh_reference_tetrahedron_is_ginh_unit_one() -> None:
    # Gmsh's reference tetrahedron is a right-corner one; it must score 1 in Gmsh's unit.
    points = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)]
    quality = tet_metrics(points, [(0, 1, 2, 3)])
    assert quality["gmsh_icn"][0] == pytest.approx(1.0, abs=1e-12)
    # ... and the regular tetrahedron is (1/sqrt2)^(2/3) = 0.7937 in that unit.
    regular = tet_metrics(regular_tet_points(), [(0, 1, 2, 3)])
    expected = (1.0 / math.sqrt(2.0)) ** (2.0 / 3.0)
    assert regular["gmsh_icn"][0] == pytest.approx(expected, abs=1e-12)
    assert expected == pytest.approx(0.7937005259840998, abs=1e-15)


# --------------------------------------------------------------------------------------
# invariance, degeneracy, orientation
# --------------------------------------------------------------------------------------
def _rotated(points: list[tuple[float, float, float]], angle: float, scale: float) -> list[tuple[float, float, float]]:
    c, s = math.cos(angle), math.sin(angle)
    out = []
    for x, y, z in points:
        x, y = scale * (c * x - s * y), scale * (s * x + c * y)
        out.append((x + 3.0, y - 2.0, z + 7.0))
    return out


def test_metrics_are_rigid_motion_and_scale_invariant() -> None:
    points = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.3, 0.9, 0.0)]
    base = tri_metrics(points, [(0, 1, 2)])
    moved = tri_metrics(_rotated(points, 0.7, 2.5), [(0, 1, 2)])
    for metric in ("min_sicn", "scaled_jacobian", "gamma", "aspect_ratio", "min_angle", "eta"):
        assert moved[metric][0] == pytest.approx(base[metric][0], abs=1e-12), metric
    assert moved["area"][0] == pytest.approx(base["area"][0] * 2.5**2, rel=1e-12)


def test_inverted_tetrahedron_is_reported_negative() -> None:
    points = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)]
    # Swap the last two vertices: same |volume|, reversed orientation.
    quality = tet_metrics(points, [(0, 1, 3, 2)])
    assert quality["volume"][0] == pytest.approx(1.0 / 6.0, rel=1e-15)
    assert quality["min_sicn"][0] < 0.0
    assert quality["gmsh_icn"][0] < 0.0
    assert quality["scaled_jacobian"][0] < 0.0


def test_flat_triangle_does_not_raise() -> None:
    points = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0)]
    quality = tri_metrics(points, [(0, 1, 2)])
    assert quality["area"][0] == 0.0
    assert quality["min_sicn"][0] == 0.0
    assert quality["gamma"][0] == 0.0
    assert quality["min_angle"][0] == pytest.approx(0.0, abs=1e-9)


def test_flat_hexahedron_is_non_positive() -> None:
    points = [
        (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
        (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
    ]
    quality = hex_metrics(points, [tuple(range(8))])
    assert quality["min_sicn"][0] <= 0.0
    assert quality["volume"][0] == pytest.approx(0.0, abs=1e-15)


# --------------------------------------------------------------------------------------
# statistics and gate
# --------------------------------------------------------------------------------------
def test_summary_stats_match_linear_interpolation() -> None:
    stats = summary_stats("min_sicn", [0.0, 1.0, 2.0, 3.0, 4.0])
    assert stats.count == 5
    assert stats.minimum == 0.0
    assert stats.median == pytest.approx(2.0)
    assert stats.p05 == pytest.approx(0.2)
    assert stats.non_positive == 1
    empty = summary_stats("min_sicn", [])
    assert math.isnan(empty.minimum) and empty.count == 0


def test_gate_passes_on_a_clean_cube_and_fails_on_a_sliver() -> None:
    cube = [
        (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 1.0), (0.0, 1.0, 1.0),
    ]
    good = measure(cube, {"hex": [tuple(range(8))]})
    assert DEFAULT_GATE.evaluate(good).passed

    sliver = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.5, 1e-4, 0.0)]
    bad = measure(sliver, {"tri": [(0, 1, 2)]})
    report = DEFAULT_GATE.evaluate(bad)
    assert not report.passed
    failing = {check.name for check in report.failures}
    assert {"tri.min_sicn", "tri.min_angle", "tri.scaled_jacobian"} <= failing
    assert "tri.edge_ratio" not in failing


def test_gate_ignores_cell_kinds_the_mesh_does_not_contain() -> None:
    points = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.5, math.sqrt(3.0) / 2.0, 0.0)]
    quality = measure(points, {"tri": [(0, 1, 2)]})
    gate = QualityGate({"hex": (MetricLimit("min_sicn", minimum=0.5),)})
    assert gate.evaluate(quality).checks == ()


def test_measure_rejects_bad_input() -> None:
    points = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)]
    with pytest.raises(ValueError):
        measure(points, {"tri": [(0, 1, 5)]})
    with pytest.raises(ValueError):
        measure(points, {"tet": [(0, 1, 2)]})
    with pytest.raises(ValueError):
        measure(points, {"wedge": [(0, 1, 0)]})
    with pytest.raises(ValueError):
        measure(points, {})


# --------------------------------------------------------------------------------------
# reference implementations (that the library is not its own judge)
# --------------------------------------------------------------------------------------
def test_triangle_reference_agrees_with_library() -> None:
    points = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.2, 0.8, 0.0), (1.0, 1.0, 0.0), (-0.4, 0.9, 0.1)]
    triangles = [(0, 1, 2), (1, 3, 2), (0, 2, 4)]
    library = tri_metrics(points, triangles)
    independent = reference.reference_triangles(points, triangles)
    for metric in independent:
        for got, ref in zip(library[metric], independent[metric]):
            if metric == "area":
                assert got == pytest.approx(ref, rel=1e-12)
            else:
                assert got == pytest.approx(ref, abs=1e-9), (metric, got, ref)


def test_tetrahedron_reference_agrees_with_library() -> None:
    points = regular_tet_points() + [(0.1, 0.2, 0.9), (1.1, 0.2, 0.9), (0.2, 1.2, 0.9), (0.4, 0.3, 1.7)]
    tets = [(0, 1, 2, 3), (4, 5, 6, 7)]
    library = tet_metrics(points, tets)
    independent = reference.reference_tets(points, tets)
    for index in range(2):
        assert library["volume"][index] == pytest.approx(independent["volume"][index], rel=1e-12)
        assert library["eta"][index] == pytest.approx(independent["eta"][index], abs=1e-10)
        # Dihedral angles: two different constructions, both exact for a regular tet.
        assert library["min_dihedral"][index] == pytest.approx(independent["min_dihedral"][index], abs=1e-9)
        assert library["max_dihedral"][index] == pytest.approx(independent["max_dihedral"][index], abs=1e-9)


def test_hexahedron_reference_finite_differences_agree_with_analytic_jacobian() -> None:
    points = [
        (0.0, 0.0, 0.0), (2.0, 0.1, 0.0), (2.2, 1.3, 0.1), (-0.1, 1.0, 0.0),
        (0.0, 0.0, 0.4), (2.0, 0.0, 0.5), (2.1, 1.1, 0.6), (0.1, 1.0, 0.4),
    ]
    library = hex_metrics(points, [tuple(range(8))])
    independent = reference.reference_hexahedra(points, [tuple(range(8))])
    assert library["scaled_jacobian"][0] == pytest.approx(independent["scaled_jacobian"][0], abs=1e-7)
    assert library["min_sicn"][0] == pytest.approx(independent["min_sicn"][0], abs=1e-7)


def test_icn_from_singular_values_matches_determinant_formula() -> None:
    pytest.importorskip("numpy")
    jac = [[2.0, 0.3, 0.1], [0.0, 1.5, 0.2], [0.0, 0.0, 0.7]]
    sj, sicn = _corner_quality(jac)
    assert reference.icn_from_singular_values(jac) == pytest.approx(abs(sicn), rel=1e-12)


def test_trimesh_is_an_independent_oracle_when_installed() -> None:
    pytest.importorskip("trimesh")
    points = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.2, 0.8, 0.0), (1.2, 1.0, 0.3)]
    triangles = [(0, 1, 2), (1, 3, 2)]
    library = tri_metrics(points, triangles)
    external = reference.trimesh_triangle_metrics(points, triangles)
    assert external is not None
    for index in range(2):
        assert library["area"][index] == pytest.approx(external["area"][index], rel=1e-12)
        assert library["min_angle"][index] == pytest.approx(external["min_angle"][index], abs=1e-9)
        assert library["max_angle"][index] == pytest.approx(external["max_angle"][index], abs=1e-9)


def test_external_availability_reports_versions() -> None:
    available = reference.external_availability()
    assert set(available) == {"gmsh", "meshmagick", "trimesh"}
    assert all(value is None or isinstance(value, str) for value in available.values())


# --------------------------------------------------------------------------------------
# exchange format and CLI
# --------------------------------------------------------------------------------------
def test_spec_round_trip(tmp_path: Path) -> None:
    spec = MeshSpec(
        points=((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)),
        elements={"tri": ((0, 1, 2),)},
        metadata={"source": "unit test"},
    )
    path = write_spec(tmp_path / "mesh.json", spec)
    loaded = load_spec(path)
    assert loaded.points == spec.points
    assert loaded.elements == spec.elements
    assert loaded.metadata == {"source": "unit test"}


def test_cli_exit_codes(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    good = write_spec(
        tmp_path / "good.json",
        MeshSpec(
            points=(
                (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
                (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 1.0), (0.0, 1.0, 1.0),
            ),
            elements={"hex": (tuple(range(8)),)},
        ),
    )
    assert cli.main(["--input", str(good), "--quiet"]) == 0
    assert cli.main(["--input", str(good), "--quiet", "--limit", "hex.min_sicn>=1.5"]) == 1
    assert cli.main(["--input", str(tmp_path / "missing.json"), "--quiet"]) == 2
    capsys.readouterr()


def test_cli_writes_json_report(tmp_path: Path) -> None:
    path = write_spec(
        tmp_path / "mesh.json",
        MeshSpec(
            points=((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.5, math.sqrt(3.0) / 2.0, 0.0)),
            elements={"tri": ((0, 1, 2),)},
        ),
    )
    target = tmp_path / "report.json"
    assert cli.main(["--input", str(path), "--json", str(target), "--quiet"]) == 0
    payload = json.loads(target.read_text(encoding="utf-8"))
    assert payload["gate"]["passed"] is True
    assert payload["quality"]["tri"]["count"] == 1
    assert payload["quality"]["tri"]["summary"]["min_sicn"]["min"] == pytest.approx(1.0, abs=1e-12)


def test_hex_jacobian_of_the_reference_cube_is_diagonal() -> None:
    cube = [
        (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 1.0), (0.0, 1.0, 1.0),
    ]
    jac = _hex_jacobian(cube, 0.0, 0.0, 0.0)
    assert jac == [[0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 0.5]]