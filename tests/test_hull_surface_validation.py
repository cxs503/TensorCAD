"""Physical validation of imported publisher hull geometry (anchor GAP-CAD-T1-IGES-PIPE).

The topology gate in :mod:`tensorcad.geomquality.trim_mesh` cannot be green for the
workshop files that actually carry hull geometry: ``KCS_hull_SVA.igs`` is a
half-breadth **open** patch set (measured: ``y in [0, B/2]``, entity census
``{128: 55}``, no cap surface at all), so ``watertight`` is false by construction
and always will be.  These tests therefore pin the *physical* claim instead -- the
geometry that comes out of the import is the published ship:

* an analytic box section reproduces its own volume/area exactly (method check);
* the KCS workshop file reproduces the published moulded displacement volume
  52 030 m3 and wetted surface 9 424 m2 within the 3 % acceptance threshold.

Both are external references; neither is a round trip through this repository's
own writer.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from tensorcad.geomquality.hull_surface import (
    PARTICULARS,
    classify_model,
    envelope_breadth,
    hull_characteristics,
    validate_against_particulars,
    _Slabs,
)

ASSET_ROOT = Path("/root/tp-agent/assets/geometry/hulls")
KCS = ASSET_ROOT / "KCS" / "KCS_hull_SVA.igs"


def _box_half_hull(length: float, half_beam: float, draft: float):
    """Half-breadth box hull, closing on ``y = 0`` and ``z = 0`` (no caps)."""
    vertices = []
    for x in (0.0, length):
        for y in (0.0, half_beam):
            for z in (0.0, draft):
                vertices.append((x, y, z))
    # vertex index = 4*xi + 2*yi + zi over (x, y, z) in ((0, L), (0, B/2), (0, z))
    faces = [
        (0, 1, 3), (0, 3, 2),      # x = 0 (transom end)
        (4, 6, 7), (4, 7, 5),      # x = L (bow end)
        (0, 1, 5), (0, 5, 4),      # y = 0 (centreplane)
        (2, 3, 7), (2, 7, 6),      # y = B/2 (side shell)
        (0, 2, 6), (0, 6, 4),      # z = 0 (bottom)
    ]
    return vertices, faces


def test_classify_half_breadth_model():
    vertices, _ = _box_half_hull(10.0, 2.5, 1.0)
    verdict = classify_model(vertices)
    assert verdict["half_breadth"] is True
    assert verdict["symmetry_factor"] == 2
    assert verdict["breadth"] == pytest.approx(2.5)


def test_classify_full_hull_model():
    vertices = [(-1.0, -3.0, 0.0), (1.0, 3.0, 0.0), (1.0, -3.0, 1.0)]
    verdict = classify_model(vertices)
    assert verdict["half_breadth"] is False
    assert verdict["symmetry_factor"] == 1


def test_envelope_breadth_is_exact_on_a_box():
    vertices, faces = _box_half_hull(10.0, 2.5, 1.0)
    slabs = _Slabs.build(vertices, faces)
    assert envelope_breadth(slabs, 5.0, 0.5) == pytest.approx(2.5)
    assert envelope_breadth(slabs, -1.0, 0.5) == 0.0  # outside the model
    assert envelope_breadth(slabs, 15.0, 0.5) == 0.0


def test_box_hull_volume_and_area_are_exact():
    length, half_beam, draft = 10.0, 2.5, 1.0
    vertices, faces = _box_half_hull(length, half_beam, draft)
    got = hull_characteristics(vertices, faces, waterline=draft, stations=40, levels=8)
    assert got.volume == pytest.approx(length * 2 * half_beam * draft, rel=1e-12)
    # side + bottom + both ends, doubled for symmetry; the centreplane adds no area
    expected_area = 2.0 * (length * draft + length * half_beam + 2.0 * half_beam * draft)
    assert got.s_wetted == pytest.approx(expected_area, rel=1e-12)
    assert got.cb == pytest.approx(1.0, rel=1e-12)
    assert got.cm == pytest.approx(1.0, rel=1e-12)
    assert got.symmetry_factor == 2


@pytest.mark.skipif(not KCS.exists(), reason=f"workshop asset absent: {KCS}")
def test_kcs_workshop_geometry_matches_published_particulars():
    """KCS IGES -> volume/Cb/S within the 3 % acceptance threshold."""
    from tensorcad.geomquality.trim_mesh import import_cad

    result = import_cad(KCS, weld_tolerance=1e-3)
    mesh = result.mesh
    particulars = PARTICULARS["KCS"]
    got = hull_characteristics(
        list(mesh.vertices), list(mesh.faces),
        waterline=particulars["draft"], stations=61, levels=48,
    )
    report = validate_against_particulars(got, particulars, tolerance_percent=3.0)

    assert classify_model(list(mesh.vertices))["half_breadth"] is True
    # The file is an open half-breadth surface model, and the report must say so
    # rather than letting a topology gate pretend the solid is closed.
    assert result.report["after_repair"]["watertight"] is False

    volume_row = report["quantities"]["volume"]
    assert abs(volume_row["error_percent"]) <= 3.0, volume_row
    assert got.beam == pytest.approx(particulars["beam"], rel=0.02)
    # Published Cb is on Lpp; the validation row puts the measured value on that basis.
    assert abs(report["quantities"]["cb"]["error_percent"]) <= 3.0, report["quantities"]["cb"]
    # `passed` is the conjunction over volume / wetted area / Cb, and it is now
    # **true**: `hull_characteristics` removes the sheets that hide under an
    # overlapping patch before summing (see `wetted_area` and
    # `test_hull_surface_dedup.py`).  The raw soup sum stays on record -- it is
    # what a naive triangle sum reports, +19.4 % on this file -- so the defect
    # remains measured rather than forgotten.
    area_row = report["quantities"]["s_wetted"]
    assert abs(area_row["error_percent"]) <= 3.0, area_row
    assert area_row["passed"] is True
    assert got.s_wetted_naive > 1.15 * particulars["s_wetted"]
    assert got.s_wetted_discarded > 0.0
    assert report["passed"] is True


@pytest.mark.skipif(not KCS.exists(), reason=f"workshop asset absent: {KCS}")
def test_kcs_station_integral_converges():
    """The sectional-area integral converges under z-grid refinement (no tuning)."""
    from tensorcad.geomquality.trim_mesh import import_cad

    result = import_cad(KCS, weld_tolerance=1e-3)
    coarse = hull_characteristics(list(result.mesh.vertices), list(result.mesh.faces),
                                  waterline=10.8, stations=41, levels=12)
    fine = hull_characteristics(list(result.mesh.vertices), list(result.mesh.faces),
                                waterline=10.8, stations=81, levels=48)
    change = abs(fine.volume - coarse.volume) / fine.volume * 100.0
    assert change < 3.0, f"volume moved {change:.2f} % between grids"

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# KVLCC2: the second of the "3 real IGES" — frame evidence is real, volume is WITHHELD
# --------------------------------------------------------------------------- #
# The KVLCC2 export (same publisher, different frame) puts z = 0 at the design
# waterline with z growing *downwards* (measured extent 0 ... -28.067 m) and x = 0
# near the bow (stern at x ~= -328 m).  Feeding it to a route that keeps
# 0 <= z <= draft keeps *nothing* — the waterline slices vanish.  That is a frame
# defect, not a shape defect, and the declarative frame transformation
# (tensorcad.fleet.frames.normalise) is the honest fix.
#
# BUT the frame fix does NOT make the volume green: mission m22903 measured the
# KVLCC2 trim-mask patch set by exact divergence theorem and got 256 657-256 983 m3
# vs the published 312 622 m3 (-17.9 %), with closure residual |∮n dA| = 11 272 m2
# (z-residual 10 003 m2 = the waterline cap stored twice).  So KVLCC2 is WITHHELD:
# the root cause is patch orientation inconsistency + waterline-cap duplication,
# not the frame.  These tests therefore pin the *frame* evidence (which is real and
# reproducible) and record the honest WITHHELD status rather than pretend a 3 % gate
# is green on a file whose exact integral is short by 17.9 %.
KVLCC2 = ASSET_ROOT / "KVLCC2" / "kvlcc2.igs"


@pytest.mark.skipif(not KVLCC2.exists(), reason=f"workshop asset absent: {KVLCC2}")
def test_kvlcc2_frame_is_z_down_stern_negative():
    """KVLCC2 IGES is a z-down, stern-negative export: normalise flips z and shifts x."""
    from tensorcad.geomquality.trim_mesh import import_cad
    from tensorcad.fleet.frames import normalise

    result = import_cad(KVLCC2, weld_tolerance=1e-3)
    mesh = result.mesh
    pts, frep = normalise(list(mesh.vertices))
    # The frame evidence is the real, reproducible claim: z-down + stern-negative.
    assert frep.z_flip is True, frep.as_dict()
    assert frep.x_shift < 0.0, frep.as_dict()   # stern-negative: translate aft-most to 0
    assert classify_model(pts)["half_breadth"] is True
    # The file is an open half-breadth surface model; the topology gate must say so.
    assert result.report["after_repair"]["watertight"] is False
