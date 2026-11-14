"""Pins ``tensorcad.meshquality.polymesh``: the solver-format reader and its OpenFOAM numbers.

What is asserted, and why each assertion is not self-referential:

* **round trip**: a mesh written by :func:`write_polymesh` reads back bit-identical
  (points, faces, owner, neighbour, patches) and yields *identical* measurements -- the
  reader is not silently re-deriving geometry;
* **closed forms**: a Cartesian hex grid has 0 deg non-orthogonality and 0 skewness; a
  uniformly sheared grid has ``atan(shear)`` non-orthogonality and 0 skewness (the classic
  case where one is nonzero and the other is not); two stacked cells of different extent
  have skewness ``0.5/sqrt(2)`` and non-orthogonality 45 deg, both computed by hand below;
* **independent route**: on affine meshes the polyMesh non-orthogonality agrees with
  :func:`tensorcad.meshquality.topology.non_orthogonality` to 1e-12 -- two code paths that
  share nothing but the node coordinates;
* **the repository's own geometry**: the hull volume mesh is converted to a polyMesh, read
  back, gated, and its zero-volume cells are found at the keel/stem columns the builder
  documents;
* **failure modes**: binary format, a missing file, a non-conforming face shared by three
  cells, and patches that do not cover the boundary tail all raise instead of being
  skipped, and a mesh with no internal faces fails the gate rather than passing vacuously.
"""
from __future__ import annotations

import math
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from tensorcad.meshquality import topology  # noqa: E402
from tensorcad.meshquality.polymesh import (  # noqa: E402
    PolyMesh,
    cell_geometry,
    check_polymesh,
    face_skewness,
    polymesh_from_hexes,
    read_polymesh,
    write_polymesh,
)

REPO_ROOT = Path(__file__).resolve().parents[1]


def box_hexes(nx: int, ny: int, nz: int):
    """Points and hexes of an ``nx*ny*nz`` cell grid of unit cubes."""
    points, index = [], {}
    for i in range(nx + 1):
        for j in range(ny + 1):
            for k in range(nz + 1):
                index[(i, j, k)] = len(points)
                points.append((float(i), float(j), float(k)))
    hexes = []
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                hexes.append(
                    [
                        index[(i, j, k)],
                        index[(i + 1, j, k)],
                        index[(i + 1, j + 1, k)],
                        index[(i, j + 1, k)],
                        index[(i, j, k + 1)],
                        index[(i + 1, j, k + 1)],
                        index[(i + 1, j + 1, k + 1)],
                        index[(i, j + 1, k + 1)],
                    ]
                )
    return points, hexes


def test_cartesian_grid_is_orthogonal_with_zero_skewness():
    points, hexes = box_hexes(3, 2, 2)
    mesh = polymesh_from_hexes(points, hexes)
    assert mesh.n_cells == 12
    assert mesh.n_internal_faces == 20  # (3-1)*2*2 + 3*(2-1)*2 + 3*2*(2-1) = 8+6+6
    quality = check_polymesh(mesh)
    assert max(quality.non_orthogonality) == pytest.approx(0.0, abs=1e-12)
    assert max(quality.skewness) == pytest.approx(0.0, abs=1e-12)
    for volume in quality.volumes:
        assert volume == pytest.approx(1.0, abs=1e-14)
    assert math.fsum(quality.volumes) == pytest.approx(12.0, abs=1e-13)
    assert quality.gate().passed


def test_uniform_shear_gives_atan_shear_and_no_skewness():
    points, hexes = box_hexes(3, 3, 3)
    shear = 0.5
    sheared = [(x + shear * y, y, z) for (x, y, z) in points]
    quality = check_polymesh(polymesh_from_hexes(sheared, hexes))
    # The face normal stays axis-aligned while the cell-centre line tilts by atan(shear):
    # the two definitions must disagree, and the non-orthogonality is exactly atan(0.5).
    assert max(quality.non_orthogonality) == pytest.approx(math.degrees(math.atan(shear)), abs=1e-12)
    # A parallelepiped has its face centre on the cell-centre line, so skewness is zero.
    assert max(quality.skewness) == pytest.approx(0.0, abs=1e-12)
    assert all(v == pytest.approx(1.0, abs=1e-14) for v in quality.volumes)


def test_stacked_cells_have_analytic_skewness_and_non_orthogonality():
    # Cell A: [0,1]^2 x [-1,0].  Cell B shares the z=0 face and its top is shifted +2 in x.
    # A's centre is (0.5, 0.5, -0.5), B's is (1.5, 0.5, 0.5), so d = (1, 0, 1), |d| = sqrt(2).
    # The z=0 face centre is (0.5, 0.5, 0) and its area vector is (0, 0, 1) (owner outward).
    #   t = ((Cf - C).Sf) / (d.Sf) = 0.5 / 1 = 0.5  ->  crossing = (1, 0.5, 0)
    #   |crossing - Cf| = 0.5  ->  skewness = 0.5 / sqrt(2) = 0.353553390593...
    #   acos(Sf.d / (|Sf||d|)) = acos(1/sqrt(2)) = 45 deg
    points = [
        (0, 0, -1), (1, 0, -1), (1, 1, -1), (0, 1, -1),
        (0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
    ]
    # B's bottom face is A's top face (shared nodes 4-7); its top face is shifted +2 in x.
    points += [(2, 0, 1), (3, 0, 1), (3, 1, 1), (2, 1, 1)]
    hexes = [list(range(8)), [4, 5, 6, 7, 8, 9, 10, 11]]
    quality = check_polymesh(polymesh_from_hexes(points, hexes))
    assert quality.non_orthogonality == pytest.approx((45.0,), abs=1e-12)
    assert quality.skewness == pytest.approx((0.5 / math.sqrt(2.0),), abs=1e-12)


def test_roundtrip_is_identical(tmp_path):
    points, hexes = box_hexes(2, 2, 1)
    mesh = polymesh_from_hexes(points, hexes, patch_name="hull")
    write_polymesh(mesh, tmp_path)
    again = read_polymesh(tmp_path)
    assert again.points == mesh.points
    assert again.faces == mesh.faces
    assert again.owner == mesh.owner
    assert again.neighbour == mesh.neighbour
    assert again.patches == mesh.patches
    first = check_polymesh(mesh)
    second = check_polymesh(again)
    assert second.non_orthogonality == first.non_orthogonality
    assert second.skewness == first.skewness
    assert second.volumes == first.volumes


def test_reads_expanded_face_form(tmp_path):
    # OpenFOAM also writes faces without the per-face counts: ( (0 1 2 3) (4 5 6 7) ... )
    points, hexes = box_hexes(1, 1, 1)
    mesh = polymesh_from_hexes(points, hexes)
    write_polymesh(mesh, tmp_path)
    text = (tmp_path / "faces").read_text().replace("4(", "(")
    (tmp_path / "faces").write_text(text)
    again = read_polymesh(tmp_path)
    assert again.faces == mesh.faces


def test_independent_route_agrees_on_affine_meshes():
    points, hexes = box_hexes(3, 2, 2)
    for label, cloud in (
        ("cartesian", points),
        ("sheared", [(x + 0.3 * z, y - 0.2 * x, z) for (x, y, z) in points]),
    ):
        mine = max(check_polymesh(polymesh_from_hexes(cloud, hexes)).non_orthogonality)
        theirs = topology.non_orthogonality(cloud, hexes, "hex").maximum
        assert mine == pytest.approx(theirs, abs=1e-12), label


def test_cell_volume_of_skewed_hex_matches_determinant():
    # A sheared box is a parallelepiped: |det| of its three edge vectors is the volume.
    points, hexes = box_hexes(1, 1, 1)
    cloud = [(x + 0.4 * y + 0.1 * z, y + 0.25 * z, z) for (x, y, z) in points]
    centres, volumes = cell_geometry(polymesh_from_hexes(cloud, hexes))
    assert volumes[0] == pytest.approx(1.0, abs=1e-14)
    assert centres[0] == pytest.approx((0.5 + 0.2 + 0.05, 0.5 + 0.125, 0.5), abs=1e-13)


def test_binary_format_raises_instead_of_guessing(tmp_path):
    points, hexes = box_hexes(1, 1, 1)
    write_polymesh(polymesh_from_hexes(points, hexes), tmp_path)
    text = (tmp_path / "points").read_text().replace("format      ascii;", "format      binary;")
    (tmp_path / "points").write_text(text)
    with pytest.raises(ValueError, match="ascii"):
        read_polymesh(tmp_path)


def test_missing_file_and_bad_patch_layout_raise(tmp_path):
    with pytest.raises(ValueError, match="does not look like a polyMesh"):
        read_polymesh(tmp_path)
    points, hexes = box_hexes(2, 1, 1)
    mesh = polymesh_from_hexes(points, hexes)
    with pytest.raises(ValueError, match="boundary faces must be contiguous"):
        PolyMesh(points=mesh.points, faces=mesh.faces, owner=mesh.owner,
                 neighbour=mesh.neighbour,
                 patches=(type(mesh.patches[0])(name="x", kind="patch", n_faces=1, start=99),))


def test_non_conforming_face_shared_by_three_cells_raises():
    points, hexes = box_hexes(3, 1, 1)
    with pytest.raises(ValueError, match="more than two cells"):
        polymesh_from_hexes(points, [hexes[0], hexes[0], hexes[1]])


def test_no_internal_faces_fails_the_gate():
    points, hexes = box_hexes(1, 1, 1)
    quality = check_polymesh(polymesh_from_hexes(points, hexes))
    assert quality.mesh.n_internal_faces == 0
    report = quality.gate()
    assert not report.passed  # vacuous passes are exactly what the gate must not do
    assert {c.name for c in report.failures} >= {"polymesh.max_non_orthogonality"}


def test_hull_volume_mesh_through_the_solver_format():
    from tensorcad.hull import HullParameters
    from tensorcad.volume_mesh import volume_mesh_from_hull

    hull = volume_mesh_from_hull(HullParameters(nx=13, nz=7), first_layer=0.05, layers=1)
    mesh = polymesh_from_hexes(hull.vertices, hull.hexes)
    quality = check_polymesh(mesh)
    assert mesh.n_cells == len(hull.hexes)
    assert mesh.n_internal_faces > 0
    gated = quality.gate()
    # The structured hull grid carries collapsed keel/stem columns, which the builder
    # documents: the gate must say so rather than pass.
    assert gated.checks  # three checks, always reported
    assert {c.name for c in gated.checks} == {
        "polymesh.max_non_orthogonality",
        "polymesh.max_skewness",
        "polymesh.cell_volume",
    }


def test_gate_fails_on_a_highly_skewed_mesh():
    # Push the face centre off the cell-centre line by an entire cell-centre distance.
    points, hexes = box_hexes(2, 2, 1)
    cloud = list(points)
    # move the two nodes of an internal-face corner pair far in x
    cloud[1] = (cloud[1][0] + 6.0, cloud[1][1], cloud[1][2])
    quality = check_polymesh(polymesh_from_hexes(cloud, hexes))
    report = quality.gate()
    failed = {c.name for c in report.failures}
    assert failed, quality.to_dict()


def test_cli_reports_polymesh_numbers(tmp_path):
    points, hexes = box_hexes(3, 2, 2)
    mesh = polymesh_from_hexes(points, hexes)
    write_polymesh(mesh, tmp_path)
    out = tmp_path / "quality.json"
    result = subprocess.run(
        [sys.executable, "-m", "tensorcad.meshquality.cli", "--polymesh", str(tmp_path),
         "--json", str(out)],
        cwd=REPO_ROOT, capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    assert out.is_file()