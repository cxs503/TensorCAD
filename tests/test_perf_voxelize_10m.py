"""10^6--10^7 triangle mesh-derived voxelisation ladder: artifact contract.

The small-scale ladder (``ladder_voxelize_mesh.json``) proved the mechanism and
carries an honest FAIL on the coarsest grid (first-order discretisation).  This
test pins the *real* 10^6--10^7 triangle ladder: the artifact must exist, carry
measured exponent gates on both size axes, keep the voxel volume within 3 % of
the analytic hull volume, and reach the 10^7 triangle band.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from tensorcad.perf_voxel_mesh import SCHEMA, check_section

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "benchmarks" / "perf_scale" / "artifacts" / "ladder_voxelize_10m.json"
DOC = ROOT / "docs" / "perf_scale_10m.md"


@pytest.fixture(scope="module")
def artifact() -> dict:
    assert ARTIFACT.exists(), (
        "ladder_voxelize_10m.json missing: run "
        "python scripts/perf_voxelize_mesh.py --bases 401:129 --levels 0,1,2,3 "
        "--grids 24:12:6,32:16:8,48:24:12,64:32:16,96:48:24 "
        f"--output {ARTIFACT}"
    )
    return json.loads(ARTIFACT.read_text(encoding="utf-8"))


def test_artifact_schema_and_verdict(artifact: dict) -> None:
    assert artifact["schema"] == SCHEMA
    assert artifact["verdict"] == "pass"
    assert artifact["gates_not_run"] == 0, "a not_run gate is not a pass"
    assert artifact["gates_checked"] >= 3


def test_reaches_10_7_triangle_band(artifact: dict) -> None:
    measured = [r for r in artifact["results"] if r["status"] == "measured"]
    assert measured, "no measured rungs"
    max_tri = max(r["triangles"] for r in measured)
    assert max_tri >= 10_000_000, (
        f"ladder tops out at {max_tri:,} triangles; the 10^7 band is the point"
    )
    # The 10^6 band must be covered too (not just a single 10^7 rung).
    assert any(1_000_000 <= r["triangles"] < 10_000_000 for r in measured)


def test_both_size_axes_have_measured_exponent_gates(artifact: dict) -> None:
    gates = {g["name"]: g for g in artifact["gates"]}
    for name in ("scaling_voxelize_mesh", "scaling_voxelize_mesh_cells"):
        assert gates[name]["status"] == "measured", (
            f"{name} is not measured: an exponent gate that cannot be fitted "
            "must be reported not_run, never pass"
        )
        assert gates[name]["passed"] is True
        assert gates[name]["ratio"] <= 1.5


def test_volume_accuracy_within_3_percent(artifact: dict) -> None:
    gate = {g["name"]: g for g in artifact["gates"]}["voxel_volume_accuracy"]
    assert gate["status"] == "measured"
    assert gate["passed"] is True
    assert gate["ratio"] <= 0.03, "voxel volume must be within 3 % of analytic"
    # Every measured rung carries its own error, none above tolerance.
    for row in artifact["results"]:
        if row["status"] == "measured" and row["volume_rel_err"] is not None:
            assert abs(row["volume_rel_err"]) <= 0.03


def test_no_failed_rungs(artifact: dict) -> None:
    failed = [r for r in artifact["results"] if r["status"] == "failed"]
    assert not failed, f"failed rungs must be fixed or documented, not committed: {failed}"


def test_doc_section_matches_artifact(artifact: dict) -> None:
    assert DOC.exists()
    problems = check_section(DOC.read_text(encoding="utf-8"), artifact)
    assert not problems, f"doc drift: {problems}"


def test_band_caps_magnitude_and_preserves_volume() -> None:
    """The narrow-band fix: magnitude capped at band, volume unchanged.

    Regression guard for the defect where ``_voxelize_rung`` computed ``band``
    and wrote it into the note but never passed it to the field -- deep-interior
    cells paid for a full nearest-surface search.
    """
    from tensorcad.hull import HullParameters, analytic_volume, generate_hull
    from tensorcad.mesh_field import MeshSDFField
    from tensorcad.voxel import to_voxel_grid

    mesh = generate_hull(HullParameters(length=10.0, beam=2.0, draft=0.8, nx=41, nz=13))
    band = 0.2
    full = MeshSDFField(mesh)
    banded = MeshSDFField(mesh, band=band)
    # a deep-interior point: the hull's z range is [-0.8, 0], so (5, 0, -0.4)
    # is mid-depth at midships, more than `band` from the bottom face.
    p = (5.0, 0.0, -0.4)
    d_full = full.signed_distance([p])[0]
    d_band = banded.signed_distance([p])[0]
    assert abs(d_band) <= band + 1e-9, "band must cap the magnitude"
    assert abs(d_full) > abs(d_band), "the full field must be deeper than the banded one"
    # the sign (parity) is unaffected by the band
    assert (d_full < 0) == (d_band < 0)
    # volume is unchanged: band >= h/2 keeps the partial-fill fraction exact
    grid = (24, 12, 6)
    v_full = to_voxel_grid(full, grid).volume()
    v_band = to_voxel_grid(banded, grid).volume()
    assert abs(v_full - v_band) < 1e-9
    # and the banded volume is still within 3 % of the analytic hull volume
    analytic = analytic_volume(10.0, 2.0, 0.8)
    assert abs(v_band - analytic) / analytic <= 0.03
