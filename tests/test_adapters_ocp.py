"""The OCP/STEP bridge, tested only where the optional wheels are installed.

Every test in this file is skipped unless ``cadquery-ocp`` imports, so the default
install (and CI) stays standard-library-plus-NumPy. That is deliberate: the STEP
bridge is the one place where a large, LGPL-licensed, GL-linked binary enters the
picture, and it must never become a condition for the geometry core's test suite.

When OCP *is* present the checks are the ones that decide whether a STEP file is worth
anything: the sewn shell has no free edges and is closed, OpenCASCADE's own checker
accepts the solid, the volume it computes matches the mesh's signed volume (which is
the parity claim -- see ``examples/step_export``), and a re-import returns the same
volume. The remaining tests pin the failure modes: a missing extra, an empty mesh, a
nonsensical scale.
"""
from __future__ import annotations

import pytest

from tensorcad.adapters.ocp import OcpUnavailable, StepExport, have_ocp, read_step, write_step
from tensorcad.geometry import Mesh
from tensorcad.geometry_field import SphereField, to_mesh
from tensorcad.mesh import total_signed_volume

requires_ocp = pytest.mark.skipif(not have_ocp(), reason="needs the optional cadquery-ocp wheels")


def _cube() -> Mesh:
    vertices = (
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (1.0, 1.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0),
        (1.0, 0.0, 1.0),
        (1.0, 1.0, 1.0),
        (0.0, 1.0, 1.0),
    )
    faces = (
        (0, 2, 1), (0, 3, 2),
        (4, 5, 6), (4, 6, 7),
        (0, 1, 5), (0, 5, 4),
        (1, 2, 6), (1, 6, 5),
        (2, 3, 7), (2, 7, 6),
        (3, 0, 4), (3, 4, 7),
    )
    return Mesh(vertices=vertices, faces=faces, units="m", coordinate_frame="unit-cube")


def test_without_the_extra_the_error_says_how_to_get_it() -> None:
    if have_ocp():
        pytest.skip("OCP is installed, so there is nothing missing to report")
    with pytest.raises(OcpUnavailable, match=r"tensorcad\[ocp\]"):
        write_step(_cube(), "unused.step")


@requires_ocp
def test_a_cube_survives_the_round_trip(tmp_path) -> None:
    path = tmp_path / "cube.step"
    record = write_step(_cube(), path, round_trip=True)
    assert isinstance(record, StepExport)
    assert record.passed
    assert record.faces == 12
    assert record.free_edges == 0
    assert record.multiple_edges == 0
    assert record.shell_closed
    assert record.solid_valid
    assert record.volume_m3 == pytest.approx(1.0, rel=1e-9)
    assert record.area_m2 == pytest.approx(6.0, rel=1e-9)
    assert path.exists() and path.stat().st_size > 0
    assert record.file_bytes == path.stat().st_size
    assert record.round_trip_volume_error == pytest.approx(0.0, abs=1e-12)


@requires_ocp
def test_the_kernel_volume_equals_the_mesh_volume(tmp_path) -> None:
    """The parity claim: sewing must not change the solid, only its representation.

    OpenCASCADE integrates the volume of the B-rep; TensorCAD integrates the divergence
    of the triangle soup. Both are exact for a closed surface, so they must agree to
    rounding -- if they do not, one of the two surfaces has a hole or a fold.
    """
    mesh = to_mesh(SphereField(0.25), resolution=(16, 16, 16), margin=2).mesh
    record = write_step(mesh, tmp_path / "sphere.step", round_trip=True)
    assert record.shell_closed and record.free_edges == 0
    assert record.volume_m3 == pytest.approx(total_signed_volume(mesh), rel=1e-12)
    assert record.round_trip_volume_error == pytest.approx(0.0, abs=1e-12)
    # An inscribed polyhedron is smaller than the analytic sphere; that error is the
    # mesher's, not the kernel's, and it is reported rather than smoothed away.
    assert record.volume_m3 < 4.0 / 3.0 * 3.141592653589793 * 0.25**3


@requires_ocp
def test_a_step_file_can_be_read_back_and_measured(tmp_path) -> None:
    path = tmp_path / "cube.step"
    write_step(_cube(), path)
    imported = read_step(path)
    assert imported.faces == 12
    assert imported.solids == 1
    assert imported.valid
    assert imported.volume_m3 == pytest.approx(1.0, rel=1e-9)


@requires_ocp
def test_two_disjoint_islands_become_two_solids(tmp_path) -> None:
    """A detached island must stay in the file: dropping it would change the volume.

    Coarse marching-tet output can leave a small closed island behind a fin root (it is
    gone by 100 cells along the hull). Sewing turns each island into its own shell, so
    the exporter writes a compound of solids and counts them, rather than making a solid
    out of the largest shell and quietly losing the rest.
    """
    cube = _cube()
    offset = (3.0, 0.0, 0.0)
    two = Mesh(
        vertices=cube.vertices + tuple((x + offset[0], y, z) for x, y, z in cube.vertices),
        faces=cube.faces + tuple(tuple(i + 8 for i in face) for face in cube.faces),
        units="m",
        coordinate_frame="unit-cube",
    )
    record = write_step(two, tmp_path / "two.step", round_trip=True)
    assert record.shells == 2
    assert record.solids == 2
    assert record.passed
    assert record.volume_m3 == pytest.approx(2.0, rel=1e-9)
    assert record.round_trip_volume_error == pytest.approx(0.0, abs=1e-12)
    assert read_step(tmp_path / "two.step").solids == 2


@requires_ocp
def test_the_result_is_json_ready(tmp_path) -> None:
    record = write_step(_cube(), tmp_path / "cube.step", round_trip=True)
    payload = record.to_dict()
    assert payload["path"].endswith("cube.step")
    assert payload["passed"] is True
    assert isinstance(payload["round_trip_volume_error"], float)


@requires_ocp
def test_the_scale_is_recorded_and_invertible(tmp_path) -> None:
    """Coordinates are written in millimetres; the record keeps the factor used."""
    path = tmp_path / "cube_mm.step"
    metric = write_step(_cube(), path, length_scale=1.0, round_trip=True)
    assert metric.length_scale == 1.0
    assert metric.volume_m3 == pytest.approx(1.0, rel=1e-9)
    assert read_step(path, length_scale=1.0).volume_m3 == pytest.approx(1.0, rel=1e-9)


@requires_ocp
def test_an_empty_mesh_is_refused(tmp_path) -> None:
    with pytest.raises(ValueError, match="empty mesh"):
        write_step(Mesh(vertices=(), faces=()), tmp_path / "empty.step")


@requires_ocp
def test_a_nonsensical_scale_is_refused(tmp_path) -> None:
    with pytest.raises(ValueError, match="length_scale"):
        write_step(_cube(), tmp_path / "cube.step", length_scale=0.0)


@requires_ocp
def test_a_nonsensical_tolerance_is_refused(tmp_path) -> None:
    with pytest.raises(ValueError, match="sew_tolerance"):
        write_step(_cube(), tmp_path / "cube.step", sew_tolerance=-1.0)


@requires_ocp
def test_reading_a_missing_file_is_an_error(tmp_path) -> None:
    with pytest.raises(FileNotFoundError):
        read_step(tmp_path / "absent.step")
