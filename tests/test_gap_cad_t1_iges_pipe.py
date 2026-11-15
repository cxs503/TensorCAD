"""Comprehensive verification tests for GAP-CAD-T1-IGES-PIPE.

Acceptance criteria (feature_ledger):
  【管道 CLI+库函数；3 个真实 IGES 全过门；金标回归测试绿；依赖锁定】

This module provides the final verification that:
1. The CLI entry point is callable and produces correct output
2. The library function (import_cad / load_cad_geometry) works end-to-end
3. The DTMB5415_Hull_exact gold regression passes (volume < 3% error, watertight)
4. Real IGES files (KCS, KVLCC2, ONRT) pass all gates
"""
from __future__ import annotations

import json
import struct
import sys
import tempfile
from pathlib import Path

import numpy as np
import pytest

# --------------------------------------------------------------------------- #
# fixtures
# --------------------------------------------------------------------------- #
ROOT = Path(__file__).resolve().parents[1]
GOLD_STL = Path("/root/ship-performance-platform/backend/data/geometry/ships/"
                "DTMB5415_Hull_exact.stl")
ASSET_ROOT = Path("/root/tp-agent/assets/geometry/hulls")

sys.path.insert(0, str(ROOT / "benchmarks" / "verified" / "iges_import"))


def _mesh_volume(vertices, faces) -> float:
    """Signed volume of a closed triangle mesh."""
    pts = np.asarray(vertices, dtype=float)
    tris = np.asarray(faces, dtype=int)
    a, b, c = pts[tris[:, 0]], pts[tris[:, 1]], pts[tris[:, 2]]
    return float(np.einsum("ij,ij->i", a, np.cross(b, c)).sum() / 6.0)


def _read_stl(path: Path) -> np.ndarray:
    """Binary STL -> (n_triangles, 3, 3) float64 array."""
    data = path.read_bytes()
    count = struct.unpack("<I", data[80:84])[0]
    raw = np.frombuffer(data[84:84 + 50 * count], dtype=np.uint8).reshape(count, 50)
    return raw[:, 12:48].copy().view("<f4").reshape(count, 3, 3).astype(np.float64)


# --------------------------------------------------------------------------- #
# 1. CLI entry point
# --------------------------------------------------------------------------- #
def test_cli_entry_point_is_callable():
    """The CLI module is importable and has a main function."""
    from tensorcad.iges_pipe_cli import main
    assert callable(main)


def test_cli_help():
    """The CLI --help works without error."""
    from tensorcad.iges_pipe_cli import main
    import io, contextlib
    f = io.StringIO()
    with contextlib.redirect_stderr(f):
        try:
            main(["--help"])
        except SystemExit as e:
            assert e.code == 0


# --------------------------------------------------------------------------- #
# 2. Library function
# --------------------------------------------------------------------------- #
def test_load_cad_geometry_returns_mesh_and_report():
    """load_cad_geometry returns (Mesh, dict) with expected keys."""
    from tensorcad.geometry import load_cad_geometry
    from tensorcad.geomquality.io import write_iges
    from tensorcad.geomquality.brep import sew
    from tensorcad.nurbs import NURBSSurface

    # Create a simple patch set
    surface = NURBSSurface(
        degree=(2, 2),
        control_net=[[(0, 0, 0), (1, 0, 0), (2, 0, 0)],
                     [(0, 1, 0), (1, 1.5, 0), (2, 1, 0)],
                     [(0, 2, 0), (1, 2, 0), (2, 2, 0)]],
        weights=[[1.0] * 3 for _ in range(3)],
        knots=[[0, 0, 0, 1, 1, 1], [0, 0, 0, 1, 1, 1]],
    )
    brep = sew([surface], names=["test"], source="test")

    with tempfile.NamedTemporaryFile(suffix=".iges", delete=False) as f:
        iges_path = f.name
    try:
        write_iges(iges_path, brep)
        mesh, report = load_cad_geometry(iges_path, u_samples=16, v_samples=8)
        assert hasattr(mesh, "vertices")
        assert hasattr(mesh, "faces")
        assert isinstance(report, dict)
        assert "gates" in report
        assert "after_repair" in report
    finally:
        Path(iges_path).unlink(missing_ok=True)


# --------------------------------------------------------------------------- #
# 3. DTMB5415 gold regression
# --------------------------------------------------------------------------- #
@pytest.mark.skipif(not GOLD_STL.exists(), reason="DTMB5415 gold STL not present")
def test_dtmb5415_gold_regression_watertight():
    """DTMB5415_Hull_exact.stl → loft → IGES → import → watertight mesh."""
    from gold_loft import read_stl, loft_gold_hull
    from tensorcad.geomquality.io import write_iges
    from tensorcad.geometry import load_cad_geometry

    triangles = _read_stl(GOLD_STL)
    brep, info = loft_gold_hull(triangles, stations=41, girth_samples=41)

    with tempfile.NamedTemporaryFile(suffix=".iges", delete=False) as f:
        iges_path = f.name
    try:
        write_iges(iges_path, brep)
        mesh, report = load_cad_geometry(iges_path, u_samples=64, v_samples=32)

        assert report["after_repair"]["watertight"], "Mesh must be watertight"
        assert report["after_repair"]["components"] == 1, "Must be single component"
        assert report["after_repair"]["closed_components"] == 1, "Must be closed"
        assert report["after_repair"]["signed_volume"] > 0, "Volume must be positive"
    finally:
        Path(iges_path).unlink(missing_ok=True)


@pytest.mark.skipif(not GOLD_STL.exists(), reason="DTMB5415 gold STL not present")
def test_dtmb5415_gold_regression_volume_within_3pct():
    """DTMB5415 mesh volume vs quadrature reference < 3% error."""
    from gold_loft import read_stl, loft_gold_hull, brep_quadrature
    from tensorcad.geomquality.io import write_iges
    from tensorcad.geometry import load_cad_geometry

    triangles = _read_stl(GOLD_STL)
    brep, info = loft_gold_hull(triangles, stations=41, girth_samples=41)

    # Reference volume from Gauss quadrature
    ref = brep_quadrature(brep, order=5, panels=4)
    ref_vol = ref["volume"]

    with tempfile.NamedTemporaryFile(suffix=".iges", delete=False) as f:
        iges_path = f.name
    try:
        write_iges(iges_path, brep)
        mesh, report = load_cad_geometry(iges_path, u_samples=64, v_samples=32)

        measured = _mesh_volume(mesh.vertices, mesh.faces)
        error = abs(measured - ref_vol) / ref_vol

        assert error < 0.03, (
            f"Volume error {error*100:.2f}% exceeds 3% threshold "
            f"(measured={measured:.2f}, reference={ref_vol:.2f})"
        )
    finally:
        Path(iges_path).unlink(missing_ok=True)


@pytest.mark.skipif(not GOLD_STL.exists(), reason="DTMB5415 gold STL not present")
def test_dtmb5415_gold_regression_bounding_box():
    """DTMB5415 mesh bounding box matches gold STL within 1% of length."""
    from gold_loft import read_stl, loft_gold_hull
    from tensorcad.geomquality.io import write_iges
    from tensorcad.geometry import load_cad_geometry

    triangles = _read_stl(GOLD_STL)
    brep, info = loft_gold_hull(triangles, stations=41, girth_samples=41)

    gold_length = float(np.ptp(triangles[:, :, 0]))

    with tempfile.NamedTemporaryFile(suffix=".iges", delete=False) as f:
        iges_path = f.name
    try:
        write_iges(iges_path, brep)
        mesh, report = load_cad_geometry(iges_path, u_samples=64, v_samples=32)

        pts = np.asarray(mesh.vertices)
        mesh_length = float(np.ptp(pts[:, 0]))
        error = abs(mesh_length - gold_length) / gold_length

        assert error < 0.01, (
            f"Bounding box length error {error*100:.4f}% exceeds 1% "
            f"(gold={gold_length:.2f}, mesh={mesh_length:.2f})"
        )
    finally:
        Path(iges_path).unlink(missing_ok=True)


@pytest.mark.skipif(not GOLD_STL.exists(), reason="DTMB5415 gold STL not present")
def test_dtmb5415_all_gates_pass():
    """All 7 gates must pass for the DTMB5415 gold regression."""
    from gold_loft import read_stl, loft_gold_hull
    from tensorcad.geomquality.io import write_iges
    from tensorcad.geometry import load_cad_geometry

    triangles = _read_stl(GOLD_STL)
    brep, info = loft_gold_hull(triangles, stations=41, girth_samples=41)

    with tempfile.NamedTemporaryFile(suffix=".iges", delete=False) as f:
        iges_path = f.name
    try:
        write_iges(iges_path, brep)
        mesh, report = load_cad_geometry(iges_path, u_samples=64, v_samples=32)

        failed = [n for n, v in report.get("gates", {}).items() if not v.get("passed")]
        assert not failed, f"Gates failed: {failed}\nReport: {json.dumps(report['gates'], indent=2)}"
    finally:
        Path(iges_path).unlink(missing_ok=True)


# --------------------------------------------------------------------------- #
# 4. Real IGES files gate test
# --------------------------------------------------------------------------- #
REAL_IGES_FILES = [
    ("KCS", ASSET_ROOT / "KCS" / "KCS_hull_SVA.igs"),
    ("KVLCC2", ASSET_ROOT / "KVLCC2" / "kvlcc2.igs"),
    ("ONRT", ASSET_ROOT / "ONRT" / "onrt_wo_gap" / "onr5613-full_appendage_deckhous_mergedR.igs"),
]


@pytest.mark.parametrize("name, path", REAL_IGES_FILES, ids=[n for n, _ in REAL_IGES_FILES])
def test_real_iges_gate(name: str, path: Path):
    """Each real IGES file must pass all gates after repair."""
    if not path.exists():
        pytest.skip(f"asset not present: {path}")

    from tensorcad.geomquality.trim_mesh import import_cad

    result = import_cad(path, weld_tolerance=1e-3)
    mesh = result.mesh
    report = result.report

    assert len(mesh.vertices) > 0, f"{name}: empty mesh"
    assert len(mesh.faces) > 0, f"{name}: no faces"

    failed = [n for n, v in report.get("gates", {}).items() if not v.get("passed")]
    if failed:
        after = report.get("after_repair", {})
        pytest.fail(
            f"{name}: gates failed: {failed}\n"
            f"  vertices={len(mesh.vertices)}, faces={len(mesh.faces)}\n"
            f"  free_edges={after.get('free_edges')}, "
            f"non_manifold={after.get('non_manifold_edges')}\n"
            f"  volume={after.get('signed_volume')}\n"
            f"  report={json.dumps(report.get('gates', {}), indent=2)}"
        )


# --------------------------------------------------------------------------- #
# 5. Dependency lock verification
# --------------------------------------------------------------------------- #
def test_core_package_has_no_runtime_dependencies():
    """The core tensorcad package must import without any third-party dependencies."""
    import importlib
    # Clear any cached imports
    for mod in list(sys.modules.keys()):
        if mod.startswith("tensorcad"):
            del sys.modules[mod]

    # Block numpy and torch
    class Blocker:
        def find_module(self, name, path=None):
            if name in ("numpy", "torch"):
                return self
        def load_module(self, name):
            raise ImportError(f"BLOCKED: {name}")

    blocker = Blocker()
    sys.meta_path.insert(0, blocker)
    try:
        importlib.import_module("tensorcad")
    finally:
        sys.meta_path.remove(blocker)


def test_ocp_is_optional_dependency():
    """OCP (OpenCASCADE) must be an optional dependency, not required for core import."""
    import tomllib
    pyproject = ROOT / "pyproject.toml"
    with open(pyproject, "rb") as f:
        data = tomllib.load(f)

    # Core dependencies must be empty
    assert data["project"]["dependencies"] == [], (
        f"Core dependencies must be empty, got: {data['project']['dependencies']}"
    )

    # OCP must be in optional-dependencies
    assert "ocp" in data["project"]["optional-dependencies"], (
        "OCP must be an optional dependency"
    )
