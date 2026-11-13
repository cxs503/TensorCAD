"""Pins the non-per-element half of ``tensorcad.meshquality``: skewness, topology, faces.

Four things are asserted, in increasing order of strength:

* the closed-form values of equiangle skewness on canonical elements -- the *definition*
  is right and its range is the documented [0, 1];
* exact agreement with the independent route in
  :mod:`tensorcad.meshquality.reference` for skewness, and with a live third-party
  implementation (``trimesh``) for every topological fact on seven surfaces, so the
  topology module is never its own judge;
* exact agreement of internal-face non-orthogonality with a closed form on a uniformly
  sheared hexahedral stack (``atan(tan(theta)) = theta``, to 1e-12), on the Cartesian
  case (0 deg) and on a sheared tetrahedral stack;
* the gate and CLI behaviour a CI job depends on: a closed hull passes, an open or
  inward-wound one fails and exit status is 1, and unreadable input is 2.

The numbers behind ``docs/mesh_quality_topology.md`` come from
``tools/mesh_quality_topology_xval.py``, which runs the same comparisons over larger
meshes and writes ``artifacts/mesh_quality_topology/results.json``.
"""
from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path

import pytest

from tensorcad.meshquality import (
    gate_non_orthogonality,
    gate_surface,
    hex_metrics,
    non_orthogonality,
    quad_metrics,
    surface_topology,
    tet_metrics,
    tri_metrics,
)
from tensorcad.meshquality.reference import (
    reference_non_orthogonality,
    reference_skewness,
    trimesh_topology,
)

# --------------------------------------------------------------------------------------
# canonical geometry
# --------------------------------------------------------------------------------------
BOX_POINTS = [
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.0),
    (1.0, 1.0, 0.0),
    (0.0, 1.0, 0.0),
    (0.0, 0.0, 1.0),
    (1.0, 0.0, 1.0),
    (1.0, 1.0, 1.0),
    (0.0, 1.0, 1.0),
]
BOX_FACES = [
    (0, 3, 2),
    (0, 2, 1),
    (4, 5, 6),
    (4, 6, 7),
    (0, 1, 5),
    (0, 5, 4),
    (1, 2, 6),
    (1, 6, 5),
    (2, 3, 7),
    (2, 7, 6),
    (3, 0, 4),
    (3, 4, 7),
]
REGULAR_TET = [
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.0),
    (0.5, math.sqrt(3.0) / 2.0, 0.0),
    (0.5, math.sqrt(3.0) / 6.0, math.sqrt(2.0 / 3.0)),
]
RIGHT_TET = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)]
UNIT_HEX = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
            (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 1.0), (0.0, 1.0, 1.0)]


def shear(x: float) -> list[list[float]]:
    """Unit hex sheared in x by ``x`` in the top face: x stays in (0, 1) for x < 1."""
    return [[px + (x if pz > 0.5 else 0.0), py, pz] for px, py, pz in UNIT_HEX]


def hex_grid(nx: int, ny: int, nz: int, shear_per_layer: float = 0.0):
    """Structured hexahedral grid, optionally sheared in x proportional to z."""

    def index(i: int, j: int, k: int) -> int:
        return k * (ny + 1) * (nx + 1) + j * (nx + 1) + i

    points = [
        (float(i) + shear_per_layer * k, float(j), float(k))
        for k in range(nz + 1)
        for j in range(ny + 1)
        for i in range(nx + 1)
    ]
    cells = [
        (
            index(i, j, k),
            index(i + 1, j, k),
            index(i + 1, j + 1, k),
            index(i, j + 1, k),
            index(i, j, k + 1),
            index(i + 1, j, k + 1),
            index(i + 1, j + 1, k + 1),
            index(i, j + 1, k + 1),
        )
        for k in range(nz)
        for j in range(ny)
        for i in range(nx)
    ]
    return points, cells


# --------------------------------------------------------------------------------------
# skewness
# --------------------------------------------------------------------------------------
def test_skewness_closed_form_triangles() -> None:
    """Equilateral 0, right isosceles 0.25, collinear 1 -- each by hand."""
    equilateral = tri_metrics(REGULAR_TET, [(0, 1, 2)])["skewness"][0]
    right = tri_metrics(RIGHT_TET, [(0, 1, 3)])["skewness"][0]
    collinear = tri_metrics([(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0)], [(0, 1, 2)])[
        "skewness"
    ][0]
    assert equilateral == pytest.approx(0.0, abs=1e-12)
    # right isosceles: 90 -> (90-60)/120 = 0.25; 45 -> (60-45)/60 = 0.25
    assert right == pytest.approx(0.25, abs=1e-12)
    assert collinear == pytest.approx(1.0, abs=1e-12)


def test_skewness_closed_form_quads_hexes_tets() -> None:
    square = quad_metrics(UNIT_HEX, [(0, 1, 2, 3)])["skewness"][0]
    cube = hex_metrics(UNIT_HEX, [(0, 1, 2, 3, 4, 5, 6, 7)])["skewness"][0]
    regular = tet_metrics(REGULAR_TET, [(0, 1, 2, 3)])["skewness"][0]
    right_corner = tet_metrics(RIGHT_TET, [(0, 1, 2, 3)])["skewness"][0]
    assert square == pytest.approx(0.0, abs=1e-12)
    assert cube == pytest.approx(0.0, abs=1e-12)
    assert regular == pytest.approx(0.0, abs=1e-12)
    # right-corner tet: three 90 deg face angles and three 45 deg ones -> max(0.25, 0.25)
    assert right_corner == pytest.approx(0.25, abs=1e-12)


def test_skewness_matches_independent_route_for_triangles() -> None:
    """The library's skewness and ``reference_skewness`` share no code path."""
    points = [
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (0.3, 0.4, 0.0),
        (0.2, 1.7, 0.5),
        (-0.4, 0.9, 1.3),
        (2.0, -0.3, 0.1),
    ]
    faces = [(0, 2, 1), (0, 1, 4), (1, 3, 4), (2, 3, 5), (0, 4, 5)]
    ours = tri_metrics(points, faces)["skewness"]
    theirs = reference_skewness(points, faces)
    for a, b in zip(ours, theirs):
        assert a == pytest.approx(b, abs=1e-12, rel=0.0)


def test_skewness_is_bounded_normalised_and_scale_invariant() -> None:
    """Equiangle skewness lives on [0, 1] (that is why it is the commerical gate metric)."""
    points = [
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (0.5, 0.02, 0.0),
        (0.5, 0.5, 0.0),
        (0.5, 0.02, 1.0),
        (2.0, 0.0, 0.0),
    ]
    faces = [(0, 2, 1), (0, 1, 3), (0, 3, 2), (1, 2, 3), (0, 1, 4), (0, 4, 5)]
    for value in tri_metrics(points, faces)["skewness"]:
        assert 0.0 - 1e-12 <= value <= 1.0 + 1e-12
    scaled = [(3.7 * x, 3.7 * y, 3.7 * z) for x, y, z in points]
    for a, b in zip(tri_metrics(points, faces)["skewness"], tri_metrics(scaled, faces)["skewness"]):
        assert a == pytest.approx(b, abs=1e-12)


def test_hex_skewness_closed_form_and_its_blind_spot() -> None:
    """Sheared hex: skewness = alpha/90 exactly. Stretched box: skewness 0 -- so it is blind
    to anisotropy, which is why ``edge_ratio`` and ``min_sicn`` are reported next to it."""
    alpha = math.degrees(math.atan(0.6))
    sheared = hex_metrics(shear(0.6), [(0, 1, 2, 3, 4, 5, 6, 7)])
    # side faces become parallelograms with angles 90 +/- alpha -> skew = alpha / 90
    assert sheared["skewness"][0] == pytest.approx(alpha / 90.0, abs=1e-12)

    stretched = [[4.0 * x, y, z] for x, y, z in UNIT_HEX]
    box = hex_metrics(stretched, [(0, 1, 2, 3, 4, 5, 6, 7)])
    assert box["skewness"][0] < 1e-12
    assert box["scaled_jacobian"][0] == pytest.approx(1.0, abs=1e-12)  # shape only
    # J = diag(4,1,1): 3 * det^(2/3) / ||J||_F^2 = 3 * 4^(2/3) / 18
    assert box["min_sicn"][0] == pytest.approx(3.0 * 4.0 ** (2.0 / 3.0) / 18.0, abs=1e-12)
    assert box["edge_ratio"][0] == pytest.approx(4.0, abs=1e-12)


# --------------------------------------------------------------------------------------
# surface topology
# --------------------------------------------------------------------------------------
def test_closed_box_topology_and_volume() -> None:
    topology = surface_topology(BOX_POINTS, BOX_FACES)
    assert topology.watertight
    assert topology.edge_manifold
    assert topology.vertex_manifold
    assert topology.winding_consistent
    assert topology.outward
    assert topology.boundary_edges == 0
    assert topology.non_manifold_edges == 0
    assert topology.duplicate_faces == 0
    assert topology.euler_characteristic == 2
    assert topology.disconnected_components == 1
    assert topology.genus == pytest.approx(0.0)
    assert topology.signed_volume == pytest.approx(1.0, abs=1e-15)
    assert topology.problems() == ()


def test_open_and_inverted_box_are_reported_not_clamped() -> None:
    open_topology = surface_topology(BOX_POINTS, BOX_FACES[:-1])
    assert not open_topology.watertight
    assert open_topology.boundary_edges == 3
    assert any("open" in problem for problem in open_topology.problems())

    inward = surface_topology(BOX_POINTS, [tuple(reversed(face)) for face in BOX_FACES])
    assert inward.watertight  # closure is unaffected by the flip
    assert inward.winding_consistent  # all faces flipped together stay consistent
    assert inward.signed_volume == pytest.approx(-1.0, abs=1e-15)
    assert not inward.outward
    assert any("point inward" in problem for problem in inward.problems())


def test_non_manifold_and_duplicate_faces_are_counted() -> None:
    # one extra triangle glued to the cube's 0-1 edge through a new apex: that edge is now
    # used by three faces and nothing else changes
    fin = surface_topology(BOX_POINTS + [(0.5, 0.5, -1.0)], BOX_FACES + [(0, 1, 8)])
    assert fin.non_manifold_edges == 1
    assert not fin.edge_manifold
    assert not fin.watertight

    duplicated = surface_topology(BOX_POINTS, BOX_FACES + [BOX_FACES[0]])
    assert duplicated.duplicate_faces == 1
    assert not duplicated.watertight


def test_two_shells_and_genus() -> None:
    points = BOX_POINTS + [(x, y, z + 3.0) for x, y, z in BOX_POINTS]
    faces = BOX_FACES + [tuple(i + 8 for i in face) for face in BOX_FACES]
    topology = surface_topology(points, faces)
    assert topology.disconnected_components == 2
    assert topology.euler_characteristic == 4
    assert topology.genus == pytest.approx(0.0)  # per shell
    assert topology.signed_volume == pytest.approx(2.0, abs=1e-15)


def test_topology_gate_reports_every_requirement() -> None:
    clean = gate_surface(surface_topology(BOX_POINTS, BOX_FACES), max_genus=0.0)
    assert [check.name for check in clean] == [
        "surface.watertight",
        "surface.edge_manifold",
        "surface.vertex_manifold",
        "surface.winding",
        "surface.outward",
        "surface.genus",
    ]
    assert all(check.passed for check in clean)

    open_checks = gate_surface(surface_topology(BOX_POINTS, BOX_FACES[:-1]), max_genus=0.0)
    failed = {check.name for check in open_checks if not check.passed}
    assert failed == {"surface.watertight", "surface.outward", "surface.genus"}


def test_topology_rejects_bad_connectivity() -> None:
    with pytest.raises(ValueError):
        surface_topology(BOX_POINTS, [(0, 1, 2, 3)])
    with pytest.raises(ValueError):
        surface_topology(BOX_POINTS, [(0, 1, 99)])


# --------------------------------------------------------------------------------------
# external oracle: trimesh
# --------------------------------------------------------------------------------------
def _oracle_cases() -> dict[str, tuple[list, list]]:
    points = list(BOX_POINTS)
    faces = list(BOX_FACES)
    two_points = BOX_POINTS + [(x, y, z + 3.0) for x, y, z in BOX_POINTS]
    two_faces = BOX_FACES + [tuple(i + 8 for i in face) for face in BOX_FACES]
    return {
        "closed_box": (points, faces),
        "open_box": (points, faces[:-1]),
        "inverted_box": (points, [tuple(reversed(face)) for face in faces]),
        "non_manifold_fin": (points, faces + [(0, 1, 6), (1, 0, 5)]),
        "duplicate_face": (points, faces + [faces[0]]),
        "two_shells": (two_points, two_faces),
        "regular_tet": (REGULAR_TET, [(0, 2, 1), (0, 1, 3), (0, 3, 2), (1, 2, 3)]),
    }


@pytest.mark.parametrize("name", sorted(_oracle_cases()))
def test_topology_agrees_with_trimesh(name: str) -> None:
    """Every topological fact must match a live third-party implementation exactly."""
    points, faces = _oracle_cases()[name]
    oracle = trimesh_topology(points, faces)
    if oracle is None:  # pragma: no cover - trimesh is installed in this environment
        pytest.skip("trimesh not importable")
    ours = surface_topology(points, faces)
    assert ours.watertight == oracle["watertight"]
    assert ours.winding_consistent == oracle["winding_consistent"]
    assert ours.euler_characteristic == oracle["euler_number"]
    assert ours.disconnected_components == oracle["body_count"]
    # trimesh's volume is signed too, so the orientation verdict survives the comparison.
    assert ours.signed_volume == pytest.approx(oracle["volume"], abs=1e-12)


def test_trimesh_agreement_covers_both_verdicts() -> None:
    """Guard against a vacuous parametrisation: the cases must disagree on watertightness."""
    verdicts = {
        surface_topology(*_oracle_cases()[name]).watertight for name in _oracle_cases()
    }
    assert verdicts == {True, False}


# --------------------------------------------------------------------------------------
# internal-face non-orthogonality
# --------------------------------------------------------------------------------------
def test_cartesian_grid_is_orthogonal() -> None:
    points, cells = hex_grid(2, 2, 2)
    faces = non_orthogonality(points, cells, "hex")
    assert faces.cells == 8
    assert faces.internal_faces == 12  # 3 internal faces per axis on a 2x2x2 grid
    assert faces.boundary_faces == 24
    assert faces.maximum == pytest.approx(0.0, abs=1e-13)
    assert gate_non_orthogonality(faces, maximum=70.0)[0].passed


@pytest.mark.parametrize("angle_deg", [0.0, 15.0, 30.0, 45.0, 60.0])
def test_sheared_hex_stack_matches_closed_form(angle_deg: float) -> None:
    """Shear x by ``tan(theta)`` per layer: the centre line tilts by exactly ``theta``."""
    points, cells = hex_grid(2, 2, 2, shear_per_layer=math.tan(math.radians(angle_deg)))
    faces = non_orthogonality(points, cells, "hex")
    assert faces.maximum == pytest.approx(angle_deg, abs=1e-12)


def test_tetrahedra_are_measured_on_the_shared_faces() -> None:
    """The six-tet Kuhn decomposition of a cube has six internal faces, all orthogonal."""
    points = [(float(x), float(y), float(z)) for x in (0, 1) for y in (0, 1) for z in (0, 1)]
    # cube corners indexed as x + 2y + 4z
    tets = [
        (0, 1, 3, 7),
        (0, 3, 2, 7),
        (0, 2, 6, 7),
        (0, 6, 4, 7),
        (0, 4, 5, 7),
        (0, 5, 1, 7),
    ]
    faces = non_orthogonality(points, tets, "tet")
    assert faces.internal_faces == 6
    assert faces.boundary_faces == 12
    # 1e-4 deg is the acos round-off floor: acos(1 - eps) ~ sqrt(2 eps) ~ 1.2e-6 deg
    assert faces.maximum == pytest.approx(0.0, abs=1e-4)


def test_non_orthogonality_gate_fails_above_the_limit() -> None:
    points, cells = hex_grid(1, 1, 3, shear_per_layer=math.tan(math.radians(75.0)))
    faces = non_orthogonality(points, cells, "hex")
    check = gate_non_orthogonality(faces, maximum=70.0)[0]
    assert not check.passed
    assert check.value == pytest.approx(75.0, abs=1e-12)
    assert faces.over(70.0)


def test_non_orthogonality_rejects_non_conforming_and_malformed_input() -> None:
    points, cells = hex_grid(1, 1, 1)
    with pytest.raises(ValueError):
        non_orthogonality(points, cells + cells + cells, "hex")  # 3 cells per face
    with pytest.raises(ValueError):
        non_orthogonality(points, [(0, 1, 2, 3, 4, 5, 6)], "hex")
    with pytest.raises(ValueError):
        non_orthogonality(points, cells, "pyramid")
    # a single cell has no internal face: reported, and the gate refuses to pass it
    single = non_orthogonality(points, cells, "hex")
    assert single.internal_faces == 0
    assert math.isnan(single.maximum)
    assert not gate_non_orthogonality(single)[0].passed


def test_duplicate_cells_are_degenerate_not_silently_orthogonal() -> None:
    """Two cells with the same nodes have zero centre offset; that must fail the gate."""
    points, cells = hex_grid(1, 1, 1)
    doubled = non_orthogonality(points, cells + cells, "hex")
    assert doubled.degenerate_faces == 6
    assert doubled.internal_faces == 0
    assert all(not math.isinf(angle) for angle in doubled.angles)
    check = gate_non_orthogonality(doubled, maximum=70.0)[0]
    assert not check.passed
    assert "degenerate" in check.detail


# --------------------------------------------------------------------------------------
# CLI: the CI contract
# --------------------------------------------------------------------------------------
def _run_cli(tmp_path: Path, payload: dict, *args: str) -> subprocess.CompletedProcess:
    spec = tmp_path / "mesh.json"
    spec.write_text(json.dumps(payload), encoding="utf-8")
    return subprocess.run(
        [sys.executable, "-m", "tensorcad.meshquality.cli", "--input", str(spec), *args],
        capture_output=True,
        text=True,
    )


def test_cli_gates_topology_and_exit_status(tmp_path: Path) -> None:
    closed = {"points": BOX_POINTS, "elements": {"tri": [list(f) for f in BOX_FACES]}}
    good = _run_cli(tmp_path, closed, "--topology", "--max-genus", "0", "--quiet")
    assert good.returncode == 0, good.stderr

    open_mesh = {"points": BOX_POINTS, "elements": {"tri": [list(f) for f in BOX_FACES[:-1]]}}
    bad = _run_cli(tmp_path, open_mesh, "--topology", "--quiet")
    assert bad.returncode == 1
    assert "check(s) failed" in bad.stderr

    volume_only = {"points": BOX_POINTS, "elements": {"hex": [[0, 1, 2, 3, 4, 5, 6, 7]]}}
    assert _run_cli(tmp_path, volume_only, "--topology").returncode == 2


def test_cli_gates_non_orthogonality_and_writes_the_report(tmp_path: Path) -> None:
    points, cells = hex_grid(1, 1, 2, shear_per_layer=math.tan(math.radians(40.0)))
    payload = {"points": points, "elements": {"hex": [list(c) for c in cells]}}
    out = tmp_path / "report.json"
    result = _run_cli(tmp_path, payload, "--max-non-orth", "30", "--json", str(out), "--quiet")
    assert result.returncode == 1
    report = json.loads(out.read_text())
    assert report["non_orthogonality"]["max_deg"] == pytest.approx(40.0, abs=1e-12)
    assert report["gate"]["failures"] == ["volume.max_non_orthogonality"]

    ok = _run_cli(tmp_path, payload, "--max-non-orth", "45", "--quiet")
    assert ok.returncode == 0


def test_non_orthogonality_matches_the_independent_arithmetic_route() -> None:
    """``acos`` of a dot product vs ``atan2`` of a cross product: must agree to 1e-12 deg."""
    points, cells = hex_grid(2, 2, 2, shear_per_layer=math.tan(math.radians(37.0)))
    ours = sorted(non_orthogonality(points, cells, "hex").angles)
    theirs = sorted(reference_non_orthogonality(points, cells, "hex"))
    assert len(ours) == len(theirs) == 12
    for a, b in zip(ours, theirs):
        assert a == pytest.approx(b, abs=1e-12)

    warped = [(x + 0.05 * math.sin(3.0 * y), y, z) for x, y, z in points]
    ours = sorted(non_orthogonality(warped, cells, "hex").angles)
    theirs = sorted(reference_non_orthogonality(warped, cells, "hex"))
    for a, b in zip(ours, theirs):
        assert a == pytest.approx(b, abs=1e-9)
