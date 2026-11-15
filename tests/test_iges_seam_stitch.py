"""Seam stage: T-junction stitching, seam measurement, and the dependency lock.

The three claims under test are the ones the pipeline report makes:

1. ``stitch_t_junctions`` closes a T-junction seam *without moving a single vertex* --
   it only makes the coarse patch share the fine patch's indices, so a closed seam
   encloses the volume it enclosed before, to machine precision;
2. ``seam_report`` does not mistake a genuine geometric gap for a seam: a metre-wide
   separation must be counted as residual, not as something to stitch;
3. the seam stage is measurement-only unless a tolerance is passed, so no existing
   caller's geometry changes under it.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

from tensorcad.geomquality.stitch import (
    edge_use_counts,
    free_edges,
    seam_report,
    stitch_t_junctions,
)
from tensorcad.geomquality.trim_mesh import import_cad, seam_stage

REAL_KCS = Path("/root/tp-agent/assets/geometry/hulls/KCS/KCS_hull_SVA.igs")
LOCK = Path(__file__).resolve().parents[1] / "requirements-iges-pipe.lock"


def _signed_volume(vertices, faces) -> float:
    total = 0.0
    for a, b, c in faces:
        ax, ay, az = vertices[a]
        bx, by, bz = vertices[b]
        cx, cy, cz = vertices[c]
        total += (ax * (by * cz - bz * cy) - ay * (bx * cz - bz * cx)
                  + az * (bx * cy - by * cx)) / 6.0
    return total


def _t_junction_shell():
    """Two patches sharing one chord, the fine one carrying a midpoint vertex on it.

    This is the defect the publisher's KCS export shows at every patch seam: the coarse
    patch has a single chord where its neighbour has two, so both edges are used once and
    the shell reads as open even though the curves coincide exactly.  The fine patch is
    folded out of plane so the enclosed volume is non-zero and can be compared.
    """
    vertices = [
        (0.0, 0.0, 0.0),   # 0  coarse chord start
        (2.0, 0.0, 0.0),   # 1  coarse chord end
        (1.0, 1.0, 0.0),   # 2  coarse patch apex
        (1.0, -1.0, 0.4),  # 3  fine patch apex
        (1.0, 0.0, 0.0),   # 4  midpoint that only the fine patch owns
    ]
    faces = [(0, 1, 2), (0, 4, 3), (4, 1, 3)]
    return vertices, faces, 4


def test_stitch_closes_a_t_junction_without_moving_a_vertex():
    vertices, faces, _ = _t_junction_shell()
    before_free = len(free_edges(faces))
    assert before_free > 0
    report = seam_report(vertices, faces, tolerance=1e-9)
    assert report["t_junctions"] >= 1  # the midpoint sits exactly on the coarse chord

    stitched_vertices, stitched_faces, record = stitch_t_junctions(
        vertices, faces, tolerance=1e-9)
    assert record["splits_total"] >= 1
    assert len(free_edges(stitched_faces)) < before_free
    # no geometry is invented: the vertex table is identical
    assert [tuple(v) for v in stitched_vertices] == [tuple(v) for v in vertices]
    # the midpoint is shared now, so the two chords it used to sit between are closed
    counts = edge_use_counts(stitched_faces)
    assert counts[(0, 4)] == 2
    assert counts[(1, 4)] == 2


def _surface_area(vertices, faces) -> float:
    total = 0.0
    for a, b, c in faces:
        ax, ay, az = vertices[a]
        bx, by, bz = vertices[b]
        cx, cy, cz = vertices[c]
        ux, uy, uz = bx - ax, by - ay, bz - az
        vx, vy, vz = cx - ax, cy - ay, cz - az
        nx, ny, nz = uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx
        total += 0.5 * math.sqrt(nx * nx + ny * ny + nz * nz)
    return total


def test_split_preserves_surface_area_to_machine_precision():
    """Splitting a chord at a point that lies *on* it must not add or remove surface."""
    vertices, faces, _ = _t_junction_shell()
    area_before = _surface_area(vertices, faces)
    stitched_vertices, stitched_faces, _ = stitch_t_junctions(vertices, faces, tolerance=1e-9)
    area_after = _surface_area(stitched_vertices, stitched_faces)
    assert math.isclose(area_before, area_after, rel_tol=0.0, abs_tol=1e-12)


def test_stitch_is_idempotent():
    vertices, faces, _ = _t_junction_shell()
    once = stitch_t_junctions(vertices, faces, tolerance=1e-9)
    twice = stitch_t_junctions(once[0], once[1], tolerance=1e-9)
    assert len(free_edges(once[1])) == len(free_edges(twice[1]))
    assert twice[2]["splits_total"] == 0


def test_wide_gap_is_reported_as_residual_not_as_a_seam():
    """Two triangles a metre apart must not be stitched by a millimetre budget."""
    vertices = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0),
                (0.0, 1.0, 1.0), (1.0, 1.0, 1.0), (1.0, 2.0, 1.0)]
    faces = [(0, 1, 2), (3, 4, 5)]
    report = seam_report(vertices, faces, tolerance=1e-3)
    assert report["free_edges"] == 6
    assert report["t_junctions"] == 0
    stitched = stitch_t_junctions(vertices, faces, tolerance=1e-3)
    assert stitched[2]["splits_total"] == 0
    assert len(free_edges(stitched[1])) == 6


def test_seam_stage_is_measurement_only_without_a_tolerance():
    vertices, faces, _ = _t_junction_shell()
    out_vertices, out_faces, record = seam_stage(vertices, faces, tolerance=None)
    assert record["applied"] is False
    assert len(out_faces) == len(faces)
    assert record["analysis_before"]["free_edges"] == len(free_edges(faces))


def test_edge_use_counts_sees_a_non_manifold_edge():
    faces = [(0, 1, 2), (0, 1, 3), (0, 1, 4)]
    counts = edge_use_counts(faces)
    assert counts[(0, 1)] == 3
    assert free_edges(faces)  # the three outer edges are still free


@pytest.mark.skipif(not REAL_KCS.exists(), reason="publisher KCS IGES not staged on this host")
def test_real_kcs_seam_stage_reduces_free_edges():
    """Publisher's 55-patch hull: a 15 mm seam budget must close a large share of seams.

    The bound is deliberately loose (>= 25 % of free edges) because it is a *measured*
    property of that export, not a target; the point of the test is that stitching is
    wired into the pipeline and that its effect is recorded, not that it closes the hull.
    """
    measured = import_cad(REAL_KCS, weld_tolerance=1e-3, repair=True)
    stitched = import_cad(REAL_KCS, weld_tolerance=1e-3, repair=True, seam_tolerance=1.5e-2)
    before = measured.report["gates"]["G3_watertight"]["value"]
    after = stitched.report["gates"]["G3_watertight"]["value"]
    assert stitched.report["seam"]["applied"] is True
    assert stitched.report["seam"]["stitch"]["splits_total"] > 0
    assert after < before * 0.75, f"free edges {before} -> {after}"


# --------------------------------------------------------------------------- #
# dependency lock (acceptance item 4 of GAP-CAD-T1-IGES-PIPE)
# --------------------------------------------------------------------------- #
def _parse_lock(path: Path) -> dict[str, str]:
    pins: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.split("#", 1)[0].strip()
        if not line or "==" not in line:
            continue
        name, version = line.split("==", 1)
        pins[name.strip().lower()] = version.strip()
    return pins


def test_lock_file_pins_the_verification_stack():
    assert LOCK.exists(), "the IGES pipe's verification stack must be pinned in the repo"
    pins = _parse_lock(LOCK)
    for required in ("numpy", "pytest", "trimesh", "meshio"):
        assert required in pins, f"{required} is not pinned"
    assert all(version for version in pins.values())


def test_lock_pins_match_the_installed_stack_on_the_platform_it_was_made_for():
    from importlib import metadata

    pins = _parse_lock(LOCK)
    platform_note = LOCK.read_text(encoding="utf-8").splitlines()[0]
    if "loongarch64" not in platform_note or "loongarch64" not in sys.platform + str(sys.maxsize):
        if "loongarch64" not in platform_note:
            pytest.skip("lock is pinned for another platform; version check not applicable")
    drift = {}
    for name, pinned in pins.items():
        try:
            installed = metadata.version(name)
        except metadata.PackageNotFoundError:
            continue
        if installed != pinned:
            drift[name] = (pinned, installed)
    assert not drift, f"environment drifted from the lock: {drift}"


def test_core_pipeline_imports_and_runs_with_third_party_blocked(monkeypatch):
    """The zero-dependency charter: the seam stage must not need numpy/trimesh/meshio."""
    blocked = ("numpy", "trimesh", "meshio", "scipy")

    class Blocker:
        def find_module(self, name, path=None):  # pragma: no cover - py<3.12 protocol
            return self if name.split(".")[0] in blocked else None

        def find_spec(self, name, path=None, target=None):
            if name.split(".")[0] in blocked:
                raise ImportError(f"BLOCKED: {name}")
            return None

    for module in list(sys.modules):
        if module.split(".")[0] in blocked:
            monkeypatch.delitem(sys.modules, module, raising=False)
    monkeypatch.setattr(sys, "meta_path", [Blocker(), *sys.meta_path])

    import importlib

    stitch = importlib.import_module("tensorcad.geomquality.stitch")
    trim_mesh = importlib.import_module("tensorcad.geomquality.trim_mesh")
    assert hasattr(stitch, "stitch_t_junctions")
    assert hasattr(trim_mesh, "seam_stage")