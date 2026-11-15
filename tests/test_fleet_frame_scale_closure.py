"""GAP-CAD-T1-BASELINE-FLEET: the two archives that were refused for *provenance*, not shape.

JBC was refused because the publisher's geometry page states no scale row, so the producer
divided the *untrimmed* breadth by B = 45 m (lambda = 39.7285, a +0.68 % length bias); ONRT
was refused with ``no_plausible_unit`` because its archive is a *unit-Lpp* model
the bare-hull sheet spans 1.00614 in x -- while the spec carried 48.9355, the workshop model
scale of a 3.147 m model, and the sampling spacing was expressed in file units, so an
"8 samples per metre" request meant eight samples over the whole hull.

These tests pin the four repairs to measurements, not to the spec file's prose:

1. ``specs_from_json`` must carry ``geometry_scale`` / ``waterline_z_file`` (it used to drop
   them, silently grading every archive as full scale);
2. the sampling spacing must be a *physical* quantity (file-unit spacing divided by
   metres-per-file-unit), so the same request means the same density on every archive;
3. JBC at the published lambda = 40 reproduces the published displacement inside the anchor's
   1 % band, and its own half-breadth x 40 reproduces B = 45 m;
4. ONRT's declared waterline is *corroborated* by the geometry (immersed depth vs published
   draft, immersed half-breadth vs B/2) -- and its residual is reported, not smoothed away:
   the hull reconstruction sits below the published displacement, so the verdict stays
   WITHHELD (an honest FAIL is not overwritten).

Run: ``python -m pytest tests/test_fleet_frame_scale_closure.py -q`` (~25 s: two archives).
"""

from __future__ import annotations

from pathlib import Path

import pytest

from tensorcad.fleet_baseline import build_ship, specs_from_json

SPECS = Path("benchmarks/verified/baseline_fleet/fleet_specs.json")
IGES = Path("artifacts/baseline_fleet/iges")
JBC = IGES / "jbc.igs"
ONRT = IGES / "onrt_wo_gap" / "onr5613-full_appendage_deckhous_mergedR.igs"

pytestmark = pytest.mark.skipif(
    not (JBC.is_file() and ONRT.is_file()),
    reason="workshop archives are not redistributed (digests only); local copy required",
)


def _spec(name: str, path: Path):
    spec = {entry.name: entry for entry in specs_from_json(SPECS)}[name]
    spec.iges_path = str(path)
    return spec


def test_loader_carries_scale_and_waterline() -> None:
    """The spec file is the single place a hull's scale lives; the loader must not drop it."""
    specs = {entry.name: entry for entry in specs_from_json(SPECS)}
    assert specs["JBC"].geometry_scale == pytest.approx(40.0)
    assert specs["KCS"].geometry_scale == pytest.approx(1.0)
    assert specs["ONRT"].geometry_scale == pytest.approx(154.0)
    assert specs["ONRT"].waterline_z_file == pytest.approx(0.0)
    assert specs["JBC"].waterline_z_file is None
    # provenance travels with the number
    assert "scale ratio of 40" in specs["JBC"].scale_basis
    assert "unit-Lpp" in specs["ONRT"].scale_basis


def test_sampling_spacing_is_physical() -> None:
    """A sample request is samples per metre of ship, not per file unit."""
    mesh, report = build_ship(_spec("JBC", JBC), dx=1.0, dz=0.25, samples_per_metre=4.0)
    source = report["source"]
    metres_per_unit = source["metres_per_file_unit_for_spacing"]
    assert metres_per_unit == pytest.approx(0.001 * 40.0, rel=1e-9)
    assert source["sampling_spacing_m"] == pytest.approx(0.25)
    assert source["sampling_spacing_file_units"] == pytest.approx(0.25 / metres_per_unit)
    assert report["mesh"]["is_watertight"] and report["mesh"]["is_edge_manifold"]


def test_jbc_verified_at_published_scale() -> None:
    """lambda = 40 closes JBC: the displacement lands inside the 1 % band."""
    mesh, report = build_ship(_spec("JBC", JBC))
    assert report["checks"]["watertight"] and report["checks"]["edge_manifold"]
    assert report["mesh"]["zero_area_faces"] == 0
    error = report["derived"]["published_volume_rel_error_percent"]
    assert abs(error) < 1.0, report["reasons"]
    assert report["verdict"] == "VERIFIED", report["reasons"]
    # the archive's own half-breadth x 40 is the published beam, which is why the factor is
    # not a fit: 562.5 mm x 40 = 45.000 m
    assert report["mesh"]["extent_m"][1] == pytest.approx(45.0, rel=0.005)
    assert report["derived"]["cb"] == pytest.approx(0.858, rel=0.02)
    # and the volume gate is graded at the publisher's own scale, not a fitted one
    assert report["source"]["scale_inferred"] is False
    assert "absolute displacement" in report["volume_gate_basis"]


def test_onrt_frame_is_declared_and_checked() -> None:
    """ONRT: the frame is declared per ship and the run carries its corroboration numbers.

    The archive is a unit-Lpp model whose superstructure sits above z = 0, so no automatic
    guard can recognise the frame; the waterline is therefore a declared data fact.  What
    the test pins is that the declaration is *in effect* (frame source, shift) and that the
    immersed depth it produces is reported next to the published draft -- the residual volume
    is then a shape/coverage fact, and the verdict stays WITHHELD instead of being promoted.
    """
    mesh, report = build_ship(_spec("ONRT", ONRT))
    frame = report["frame_check"]
    assert frame is not None
    assert frame["frame_source"] == "declared_waterline_z_file"
    assert frame["waterline_z_file"] == 0.0
    assert frame["frame_shift_z_m"] == pytest.approx(report["spec"]["draft"])
    # the two measured corroboration numbers travel with the verdict
    assert frame["immersed_depth_m"] > 0.0
    assert isinstance(frame["draft_relative_error_percent"], float)
    assert isinstance(frame["beam_relative_error_percent"], float)
    assert report["checks"]["watertight"] and report["checks"]["edge_manifold"]
    # honest FAIL: the residual is reported and the verdict is not promoted
    assert report["verdict"] == "WITHHELD"
    assert any("volume" in reason for reason in report["reasons"])


def test_onrt_wrong_waterline_is_refused() -> None:
    """A declared waterline that the geometry does not support is a refusal, not a grade."""
    spec = _spec("ONRT", ONRT)
    spec.waterline_z_file = -0.02  # 3 m below the real waterline: cuts off the topsides
    mesh, report = build_ship(spec, dx=1.0, dz=0.25)
    assert report["frame_check"]["ok"] is False
    assert report["verdict"] == "WITHHELD"
    assert any("not corroborated" in reason for reason in report["reasons"])