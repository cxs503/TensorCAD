"""Tests for the baseline fleet gate: placeholders must fail, real hulls must pass.

The fixtures are synthetic (a triangulated unit cube, a coarse box "hull" and a
deliberately duplicated pair), so the suite is offline and independent of the
workshop downloads in ``artifacts/baseline_fleet``.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from tensorcad.fleet_gate import (
    MIN_HULL_TRIANGLES,
    PLACEHOLDER,
    VERIFIED,
    WITHHELD,
    audit_fleet,
    audit_mesh,
    judge_mesh,
    mesh_digest,
    read_stl,
    validate_source,
)
from tensorcad.geometry import Mesh
from tensorcad.mesh import is_edge_manifold, is_watertight, total_signed_volume

CUBE_VERTICES = [
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.0),
    (1.0, 1.0, 0.0),
    (0.0, 1.0, 0.0),
    (0.0, 0.0, 1.0),
    (1.0, 0.0, 1.0),
    (1.0, 1.0, 1.0),
    (0.0, 1.0, 1.0),
]
CUBE_FACES = [
    (0, 2, 1), (0, 3, 2),          # bottom, outward -z
    (4, 5, 6), (4, 6, 7),          # top
    (0, 1, 5), (0, 5, 4),          # -y
    (1, 2, 6), (1, 6, 5),          # +x
    (2, 3, 7), (2, 7, 6),          # +y
    (3, 0, 4), (3, 4, 7),          # -x
]


def _write_binary_stl(path: Path, mesh: Mesh) -> Path:
    import struct

    with path.open("wb") as handle:
        handle.write(b"tensorcad fleet_gate test".ljust(80, b" "))
        handle.write(struct.pack("<I", len(mesh.faces)))
        for face in mesh.faces:
            handle.write(struct.pack("<3f", 0.0, 0.0, 0.0))
            for index in face:
                handle.write(struct.pack("<3f", *mesh.vertices[index]))
            handle.write(struct.pack("<H", 0))
    return path


def _subdivided_cube(levels: int) -> Mesh:
    """Cube with ``4**levels`` triangles per face -- enough to clear the facet gate."""
    grid = 2**levels + 1
    vertices: list[tuple[float, float, float]] = []
    faces: list[tuple[int, int, int]] = []
    axis = [(0, 1, 2), (1, 2, 0), (2, 0, 1)]
    for sign in (0.0, 1.0):
        for base in axis:
            origin = len(vertices)
            for i in range(grid):
                for j in range(grid):
                    point = [0.0, 0.0, 0.0]
                    point[base[0]] = sign
                    point[base[1]] = i / (grid - 1)
                    point[base[2]] = j / (grid - 1)
                    vertices.append((point[0], point[1], point[2]))
            flip = bool(sign)
            for i in range(grid - 1):
                for j in range(grid - 1):
                    a = origin + i * grid + j
                    b, c, d = a + 1, a + grid, a + grid + 1
                    if flip:
                        faces += [(a, d, b), (a, c, d)]
                    else:
                        faces += [(a, b, d), (a, d, c)]
    return Mesh(vertices=vertices, faces=faces)


def test_mesh_helpers_agree_on_synthetic_cube() -> None:
    mesh = Mesh(vertices=CUBE_VERTICES, faces=CUBE_FACES)
    assert is_watertight(mesh) and is_edge_manifold(mesh)
    assert math.isclose(abs(total_signed_volume(mesh)), 1.0, rel_tol=1e-12)


def test_binary_and_ascii_stl_round_trip(tmp_path: Path) -> None:
    mesh = Mesh(vertices=CUBE_VERTICES, faces=CUBE_FACES)
    binary = _write_binary_stl(tmp_path / "cube.stl", mesh)
    assert len(read_stl(binary).faces) == len(mesh.faces)
    ascii_path = tmp_path / "cube_ascii.stl"
    lines = ["solid cube"]
    for face in mesh.faces:
        lines.append("  facet normal 0 0 0")
        lines.append("    outer loop")
        for index in face:
            x, y, z = mesh.vertices[index]
            lines.append(f"      vertex {x!r} {y!r} {z!r}")
        lines.append("    endloop")
        lines.append("  endfacet")
    lines.append("endsolid cube")
    ascii_path.write_text("\n".join(lines) + "\n")
    assert len(read_stl(ascii_path).faces) == len(mesh.faces)


def test_placeholder_box_is_rejected_by_facet_count(tmp_path: Path) -> None:
    """304 facets cannot carry a hull form -- the gate must say so, not pass it."""
    coarse = Mesh(
        vertices=[(0.0, 0.0, 0.0), (20.0, 0.0, 0.0), (20.0, 4.0, 0.0), (0.0, 4.0, 0.0),
                  (0.0, 0.0, 2.0), (20.0, 0.0, 2.0), (20.0, 4.0, 2.0), (0.0, 4.0, 2.0)],
        faces=CUBE_FACES,
    )
    audit = audit_mesh(_write_binary_stl(tmp_path / "KCS_Hull.stl", coarse))
    assert audit.triangles < MIN_HULL_TRIANGLES
    verdict = judge_mesh("KCS", audit, published_volume=52_030.0)
    assert verdict.grade == PLACEHOLDER
    assert "triangles" in verdict.reasons[0]
    assert verdict.volume_error is None


def test_duplicate_bytes_across_two_ships_is_placeholder(tmp_path: Path) -> None:
    mesh = _subdivided_cube(4)
    path = _write_binary_stl(tmp_path / "KVLCC2_Hull.stl", mesh)
    twin = _write_binary_stl(tmp_path / "KVLCC1_Hull.stl", mesh)
    assert mesh_digest(read_stl(path)) == mesh_digest(read_stl(twin))
    report = audit_fleet(tmp_path, published_volumes={"KVLCC2": 312_622.0, "KVLCC1": 312_622.0})
    grades = {ship["ship"]: ship for ship in report["ships"]}
    assert grades["KVLCC1"]["grade"] == PLACEHOLDER
    assert grades["KVLCC2"]["grade"] == PLACEHOLDER
    assert "byte-identical" in " ".join(grades["KVLCC2"]["reasons"])
    assert report["pass"] is False


def test_watertight_form_passes_volume_gate_and_open_form_withholds(tmp_path: Path) -> None:
    solid = _subdivided_cube(4)
    closed_path = _write_binary_stl(tmp_path / "Good_Hull.stl", solid)
    audit = audit_mesh(closed_path)
    assert audit.watertight and audit.components == 1
    ok = judge_mesh("GOOD", audit, published_volume=abs(audit.volume) * 1.005)
    assert ok.grade == VERIFIED and ok.volume_error is not None and ok.volume_error < 0.01
    bad = judge_mesh("GOOD", audit, published_volume=abs(audit.volume) * 1.05)
    assert bad.grade == WITHHELD and bad.volume_error > 0.01

    open_mesh = Mesh(vertices=solid.vertices, faces=solid.faces[:-2])
    open_audit = audit_mesh(_write_binary_stl(tmp_path / "Open_Hull.stl", open_mesh))
    assert not open_audit.watertight
    assert judge_mesh("OPEN", open_audit, published_volume=1.0).grade == WITHHELD


def test_provenance_validator_requires_official_host_licence_and_sha() -> None:
    good = {
        "source_url": "https://simman2020.kr/download/KVLCC2-Hull.zip",
        "publisher": "SIMMAN 2020 Workshop",
        "license": "workshop benchmark data, research use",
        "sha256": "dc032647f00a68d3" + "0" * 48,
        "retrieved_at": "2026-11-13T02:56:00+00:00",
    }
    assert validate_source(good) == []
    assert "missing license" in validate_source({**good, "license": ""})
    assert "missing sha256" in validate_source({**good, "sha256": ""})
    assert any("not a known workshop publisher" in p
               for p in validate_source({**good, "source_url": "https://example.com/kvlcc2.zip"}))
    assert any("64 lowercase hex" in p for p in validate_source({**good, "sha256": "abc"}))


def test_audit_fleet_writes_json_report(tmp_path: Path) -> None:
    mesh = _subdivided_cube(4)
    _write_binary_stl(tmp_path / "A_Hull.stl", mesh)
    report = audit_fleet(tmp_path, published_volumes={"A": abs(total_signed_volume(mesh))})
    assert report["schema"].endswith("/1.0")
    assert json.dumps(report)  # serializable
    assert report["counts"].get(VERIFIED) == 1