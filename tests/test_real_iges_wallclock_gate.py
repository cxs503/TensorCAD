"""Wall-clock + gate regression for the three real publisher IGES files.

Anchor GAP-CAD-T1-IGES-PIPE, acceptance criterion 2 ("3 real IGES pass all gates").
Rounds 1-6 left this measurement open (m23992 §2.2): after the 142 field-position
fix and the NURBS derivative memoisation, nobody had yet timed KVLCC2/ONRT on a
node with asset access.  m24626 closed the measurement (see
``docs/gap_cad_t1_iges_pipe_wallclock_m24626.md``) and this test pins the result:

* each file imports through the trimmed route within the 120 s per-file budget;
* G1 (all faces read), G2 (trim on surface) and G5 (positive volume) pass.

The full watertight verdict (G3/G4/G6) is *not* claimed for publisher files: the
trimmed route tessellates each face independently and the inter-face seams are not
stitched; the watertight hulls come from the fleet_baseline half-breadth route.
Files absent from the asset store are skipped, matching the existing
``test_real_iges_gate.py`` convention.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "assets" / "geometry_registry.json"

# Per-file test budget established by the earlier rounds (120 s per file).
WALLCLOCK_BUDGET_SECONDS = 120.0

# The three real publisher files acceptance criterion 2 names.
SPECIMENS = {
    "KCS": "KCS_hull_SVA.igs",
    "KVLCC2": "kvlcc2.igs",
    "ONRT": "onr5613-full_appendage_deckhous_mergedR.igs",
}


def _asset(name: str) -> Path | None:
    registry = json.loads(REGISTRY.read_text())
    store = Path(registry.get("store_root", ""))
    for asset in registry["assets"]:
        if asset["asset_id"].endswith(name):
            path = store / asset["uri"][len("asset://"):]
            return path if path.exists() else None
    return None


@pytest.mark.parametrize("hull", sorted(SPECIMENS))
def test_real_iges_within_wallclock_budget_and_core_gates(hull: str) -> None:
    source = _asset(SPECIMENS[hull])
    if source is None:
        pytest.skip("publisher geometry not present in the asset store")
    from tensorcad.geomquality.trim_mesh import import_cad

    started = time.time()
    result = import_cad(source)
    elapsed = time.time() - started
    report = result.report
    gates = report.get("gates", {})

    # The budget is the point of this test: a regression that pushes the import
    # back over 120 s must fail loudly, not time out the whole suite.
    assert elapsed <= WALLCLOCK_BUDGET_SECONDS, (
        f"{hull}: import took {elapsed:.1f} s > {WALLCLOCK_BUDGET_SECONDS:.0f} s budget"
    )
    # G1: every face in the file was read and tessellated (none refused).
    assert gates["G1_all_faces_read"]["passed"], report["import"].get("refused")
    # G2: every trim loop inverts onto its surface within the transfer budget.
    assert gates["G2_trim_on_surface"]["passed"], gates["G2_trim_on_surface"]
    # G5: the tessellated model has positive signed volume.
    assert gates["G5_positive_volume"]["passed"], gates["G5_positive_volume"]
