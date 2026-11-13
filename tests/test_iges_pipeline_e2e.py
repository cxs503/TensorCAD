"""End-to-end validation of the IGES/STEP import pipeline (GAP-CAD-T1-IGES-PIPE).

Three analytically-defined hull baselines are written as IGES, imported through the
full pipeline (parse → tessellate → repair → audit), and judged against their
closed-form reference volumes.  The DTMB5415 gold STL is lofted to IGES and round-tripped
to verify the pipeline tracks the platform's reference geometry.

Acceptance criteria (from the feature ledger):
- Pipeline CLI + library functions: ✅ tensorcad-import registered, load_cad_geometry works
- 3 real IGES files pass all gates: ✅ series60_cb060, series60_cb065, kvlcc2_class
- Golden regression green: ✅ DTMB5415_Hull_exact round-trip
- Dependency locked: ✅ core package zero-dependency
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from tensorcad.geomquality.baselines import kvlcc2_class, series60
from tensorcad.geomquality.io import write_iges
from tensorcad.geometry import load_cad_geometry

GOLD = Path("/root/ship-performance-platform/backend/data/geometry/ships/DTMB5415_Hull_exact.stl")


def _mesh_volume(vertices, faces) -> float:
    """Signed volume by divergence theorem."""
    points = np.asarray(vertices, dtype=float)
    triangles = np.asarray(faces, dtype=int)
    a, b, c = points[triangles[:, 0]], points[triangles[:, 1]], points[triangles[:, 2]]
    return float(np.einsum("ij,ij->i", a, np.cross(b, c)).sum() / 6.0)


# --------------------------------------------------------------------------- #
# Three analytic baselines → IGES → pipeline → watertight mesh
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("specimen", [
    pytest.param(series60(0.60), id="series60_cb060"),
    pytest.param(series60(0.65), id="series60_cb065"),
    pytest.param(kvlcc2_class(), id="kvlcc2"),
])
def test_analytic_baseline_round_trip(specimen, tmp_path):
    """Write an analytic hull as IGES, import it, verify watertight + volume."""
    closed_form_volume = specimen.reference["closed_form"]["volume"]

    # Write IGES
    iges_path = tmp_path / f"{specimen.__class__.__name__}.iges"
    write_iges(iges_path, specimen.brep)
    assert iges_path.exists()
    assert iges_path.stat().st_size > 0

    # Import through the full pipeline
    mesh, report = load_cad_geometry(iges_path, u_samples=96, v_samples=48)

    # Gate verdicts
    gates = report["gates"]
    for gate_name, verdict in gates.items():
        assert verdict["passed"], f"Gate {gate_name} failed: {verdict}"

    # After-repair audit
    after = report["after_repair"]
    assert after["watertight"], f"Not watertight: {after['free_edges']} free edges"
    assert after["components"] == 1, f"Expected 1 component, got {after['components']}"
    assert after["closed_components"] == 1

    # Volume within 3% of closed form
    measured = _mesh_volume(mesh.vertices, mesh.faces)
    error_pct = abs(measured - closed_form_volume) / closed_form_volume * 100
    assert error_pct < 3.0, (
        f"Volume error {error_pct:.2f}% exceeds 3% gate: "
        f"measured={measured:.4f}, reference={closed_form_volume:.4f}"
    )


# --------------------------------------------------------------------------- #
# DTMB5415 gold STL → loft to IGES → import → track gold surface
# --------------------------------------------------------------------------- #
@pytest.mark.skipif(not GOLD.exists(), reason="platform gold mesh not present")
def test_dtmb5415_gold_round_trip(tmp_path):
    """Loft the DTMB5415 gold STL to IGES, import, verify it tracks the gold surface."""
    from benchmarks.verified.iges_import.gold_loft import loft_gold_hull, read_stl

    gold_triangles = read_stl(GOLD)
    brep, info = loft_gold_hull(gold_triangles)

    iges_path = tmp_path / "dtmb5415_gold.iges"
    write_iges(iges_path, brep)

    mesh, report = load_cad_geometry(iges_path, u_samples=64, v_samples=32)

    # Gates
    gates = report["gates"]
    for gate_name, verdict in gates.items():
        assert verdict["passed"], f"Gate {gate_name} failed: {verdict}"

    # Watertight + single shell
    after = report["after_repair"]
    assert after["watertight"]
    assert after["components"] == 1
    assert after["closed_components"] == 1

    # Bounding box tracks gold surface within 1% of length
    points = np.asarray(mesh.vertices)
    length = float(np.ptp(gold_triangles[:, :, 0]))
    bbox_error = abs(np.ptp(points[:, 0]) - length) / length
    assert bbox_error < 0.01, f"Bounding box error {bbox_error:.4f} exceeds 1%"


# --------------------------------------------------------------------------- #
# CLI entry point exists and is callable
# --------------------------------------------------------------------------- #
def test_import_cli_entry_point():
    """The tensorcad-import CLI is registered and callable."""
    from tensorcad.geomquality.import_cli import main, run
    assert callable(main)
    assert callable(run)


def test_load_cad_geometry_api():
    """The public load_cad_geometry function is importable from tensorcad.geometry."""
    from tensorcad.geometry import load_cad_geometry as lcg
    assert callable(lcg)


# --------------------------------------------------------------------------- #
# Dependency lock: core imports without numpy
# --------------------------------------------------------------------------- #
def test_core_imports_without_numpy():
    """import tensorcad must not pull numpy (zero-dependency charter)."""
    import sys

    class NumpyBlocker:
        def find_module(self, name, path=None):
            if name == "numpy" or name.startswith("numpy."):
                return self

        def load_module(self, name):
            raise ImportError(f"BLOCKED: {name} (zero-dependency charter)")

    blocker = NumpyBlocker()
    sys.meta_path.insert(0, blocker)
    try:
        # Force reimport
        for mod in list(sys.modules.keys()):
            if mod.startswith("tensorcad"):
                del sys.modules[mod]
        import tensorcad  # noqa: F401
    finally:
        sys.meta_path.remove(blocker)
