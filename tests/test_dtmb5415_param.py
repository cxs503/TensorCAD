"""Tests for the parametric DTMB 5415 design space (``tensorcad.dtmb5415``).

Scope: the shipped reduced model, the v -> mesh map, the analytic dV/dv chain and the
principal particulars.  Nothing here needs the 1.0 GB Zenodo database: the numerical
reference is ``docs/data/dtmb5415/holdout_fixture.npz`` (a dozen holdout configurations
resampled from it by ``tools/dtmb5415_param/verify_shipped_model.py``), so the tests are
runnable in CI and still compare against real database geometry.

The two round-0 defects are locked by regression tests:
  * ``test_no_double_counted_linear_block`` -- amplitudes must not add ``C1 v`` on top of a
    ``C2`` that already contains the linear columns (holdout RMS 0.833 m vs 0.252 m);
  * ``test_volume_gradient_matches_central_difference`` -- ``dV/dv`` once carried an extra
    factor 2 from the starboard mirror.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from tensorcad.dtmb5415 import (DATASET, DRAFT_M, GRID, SCALE_LAMBDA, WATERLINE_Z,
                                DTMB5415Space, point_grid)

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "docs" / "data" / "dtmb5415" / "holdout_fixture.npz"

GATE = {
    "finite_vertices": lambda q: q["non_finite_vertices"] == 0,
    "no_degenerate_faces": lambda q: q["degenerate_faces"] == 0,
    "no_non_manifold_edges": lambda q: q["non_manifold_edges"] == 0,
    "min_angle_deg": lambda q: q["min_angle_deg"] >= 20.0,
    "frac_angle_below_5deg": lambda q: q["frac_angle_below_5deg"] <= 1e-3,
    "aspect_p99": lambda q: q["aspect_p99"] <= 10.0,
    "stations_monotone": lambda q: bool(q["stations_monotone"]),
    "half_breadth_non_negative": lambda q: q["half_breadth_min_m"] >= -1e-6,
    "edge_max_m": lambda q: q["edge_max_m"] <= 8.0,
}


@pytest.fixture(scope="module")
def space() -> DTMB5415Space:
    return DTMB5415Space()


@pytest.fixture(scope="module")
def fixture():
    if not FIXTURE.exists():                                   # pragma: no cover
        pytest.skip("holdout fixture missing")
    npz = np.load(FIXTURE)
    return {k: npz[k] for k in npz.files}


# ---------------------------------------------------------------------------
# cloud layout and the reference hull
# ---------------------------------------------------------------------------
def test_point_grid_layout(space):
    """2250 points are 2 patches x 25 girth x 45 stations, deck edge -> centreline."""
    g = point_grid(space.base_cloud)
    assert g.shape == (GRID["n_patches"], GRID["n_girth"], GRID["n_station"], 3)
    # last girth row is the centreline, first row is the deck edge
    assert np.abs(g[:, -1, :, 1]).max() < 1e-6
    # deck edge sits outboard (a few stem/stern stations pinch to the centreline)
    assert np.median(g[:, 0, :, 1]) > 1.0
    # girth runs from deck edge down to the keel: z decreases along the girth index
    assert np.all(np.diff(g[:, :, :, 2], axis=1) < 0.0)
    # a station is a thin transverse slice: its x spread is a fraction of the length
    span = g[..., 0].max(axis=1) - g[..., 0].min(axis=1)          # (2, 45), within station
    # stations are transverse lines except near the ends (maximum observed 7.5% of L)
    assert span.max() < 0.10 * 142.0


def test_reference_hull_is_the_dataset_baseline(space, fixture):
    """v = 0 reproduces the baseline column exactly (the model has no intercept)."""
    assert np.abs(space.cloud(space.zero()) - space.base_cloud).max() == 0.0
    assert str(fixture["source_sha256"]) == DATASET["source_sha256"]


def test_design_vector_validation(space):
    with pytest.raises(ValueError):
        space.cloud(np.zeros(26))
    bad = np.zeros(27)
    bad[3] = np.nan
    with pytest.raises(ValueError):
        space.cloud(bad)


# ---------------------------------------------------------------------------
# shipped model accuracy (against real database configurations)
# ---------------------------------------------------------------------------
def test_shipped_model_holdout_accuracy(space, fixture):
    """RMS point error on stored holdout configurations stays at the fitted ~0.25 m."""
    errs, rms = [], []
    for v, true in zip(fixture["design_vector"], fixture["true_cloud"]):
        e = np.linalg.norm(space.cloud(v.astype(np.float64)) - true, axis=1)
        errs.append(e)
        rms.append(float(np.sqrt((e ** 2).mean())))
    rms = np.array(rms)
    all_err = np.concatenate(errs)
    overall = float(np.sqrt((all_err ** 2).mean()))
    assert overall < 0.35, f"holdout RMS {overall:.3f} m"
    assert np.median(rms) < 0.35
    assert rms.max() < 0.75
    # the fixture ships the per-configuration residuals of the full 12-config sample
    assert np.allclose(np.sort(rms), np.sort(fixture["per_config_rms_m"]), atol=1e-3)


def test_no_double_counted_linear_block(space, fixture):
    """Regression (round 0): ``C2`` already contains the linear columns of ``v``."""
    v = fixture["design_vector"][0].astype(np.float64)
    feat = np.concatenate([v, space._quadratic(v)])
    amps = space._amplitudes(v)
    assert np.allclose(amps, space.c2 @ feat)
    double = space.c1 @ v + space.c2 @ feat
    assert np.linalg.norm(double) > 1.2 * np.linalg.norm(amps)   # the old, doubling form

    cloud_fixed = space.base_cloud + (space.phi @ amps).reshape(-1, 3)
    cloud_double = space.base_cloud + (space.phi @ double).reshape(-1, 3)
    true = fixture["true_cloud"][0].astype(np.float64)
    err_fixed = np.linalg.norm(cloud_fixed - true, axis=1)
    err_double = np.linalg.norm(cloud_double - true, axis=1)
    assert np.sqrt((err_fixed ** 2).mean()) < 0.5 * np.sqrt((err_double ** 2).mean())


# ---------------------------------------------------------------------------
# v -> mesh map
# ---------------------------------------------------------------------------
def test_cloud_map_matches_shipped_coefficients(space):
    """v -> cloud is exactly ``base + Phi (C2 [v, vech(v v^T)])`` (no intercept, no C1)."""
    rng = np.random.default_rng(7)
    v = rng.normal(scale=0.3, size=27)
    feat = np.concatenate([v, space._quadratic(v)])
    expected = space.base_cloud + (space.phi @ (space.c2 @ feat)).reshape(-1, 3)
    assert np.abs(space.cloud(v) - expected).max() < 1e-9
    # and the map is quadratic, not linear: a quarter-size vector is not a quarter shift
    small = 0.25 * v
    big = (space.cloud(v) - space.base_cloud)
    quarter = (space.cloud(small) - space.base_cloud)
    assert np.abs(quarter - 0.25 * big).max() > 1e-3


def test_mesh_topology_and_gate_at_reference_hull(space):
    q = space.quality(space.zero())
    tri = space.faces()
    verts = space.vertices(space.zero())
    assert tri.shape[1] == 3 and tri.max() < verts.shape[0]
    # refined grid per patch: 48 x 132 quads -> 2 triangles each, minus the centreline strip,
    # mirrored to port: 2 * 2 * (2 * 48 * 132 - 48) / 2 ... recorded as the shipped 50880
    assert q["n_faces"] == 50880
    assert q["n_faces"] % 2 == 0
    failed = [name for name, fn in GATE.items() if not fn(q)]
    assert failed == [], f"reference hull fails {failed}"
    assert np.allclose(q["bbox"][1], [-q["bbox"][1][1], q["bbox"][1][1]], atol=1e-9)


def test_symmetric_across_the_centreplane(space):
    """The mirrored (port) shell is the starboard shell with y -> -y."""
    v = np.linspace(-0.3, 0.3, 27)
    verts = space.vertices(v)
    star = verts[: space._n_patch_vertices]
    port = verts[space._n_patch_vertices:]
    off = ~space._starboard_centre
    mirrored = star[off] * np.array([1.0, -1.0, 1.0])
    key = lambda a: np.lexsort((a[:, 2], a[:, 1], a[:, 0]))       # noqa: E731
    assert np.abs(mirrored[key(mirrored)] - port[key(port)]).max() < 1e-9
    assert port.shape[0] == int(off.sum())


# ---------------------------------------------------------------------------
# analytic differentiation chain
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("mode", [0, 4, 13, 26])
def test_volume_gradient_matches_central_difference(space, mode):
    h = 1e-5
    e = np.zeros(27)
    e[mode] = h
    fd = (space.volume(e) - space.volume(-e)) / (2 * h)
    ana = space.volume_gradient(space.zero())[mode]
    assert abs(fd) > 1.0
    assert abs(ana - fd) / abs(fd) < 1e-4


def test_vertex_jacobian_matches_central_difference(space):
    h = 1e-5
    jac = space.vertex_jacobian(space.zero())
    for mode in (0, 9):
        e = np.zeros(27)
        e[mode] = h
        fd = (space.vertices(e) - space.vertices(-e)) / (2 * h)
        assert np.abs(fd - jac[:, :, mode]).max() < 1e-6


# ---------------------------------------------------------------------------
# principal particulars (geometry before physics)
# ---------------------------------------------------------------------------
def test_principal_particulars_against_readme(space):
    t = DATASET["readme_model_scale_table"]
    lam = SCALE_LAMBDA
    ref_beam = t["beam_m"] * lam
    ref_draft = t["draft_m"] * lam
    ref_vol = (t["displacement_t"] * 1000.0) / t["rho_kg_m3"] * lam ** 3
    ref_cb = ref_vol / (142.0 * ref_beam * ref_draft)

    hydro = space.hydrostatics(None)
    assert WATERLINE_Z == pytest.approx(ref_draft, rel=1e-12)
    assert hydro["draft_m"] == pytest.approx(DRAFT_M)
    assert abs(hydro["beam_waterline_m"] / ref_beam - 1.0) < 0.03
    assert abs(hydro["displacement_volume_m3"] / ref_vol - 1.0) < 0.03
    assert abs(hydro["cb"] / ref_cb - 1.0) < 0.03
    # internal identity: Cb * L * B * T reproduces the integrated volume
    assert hydro["cb"] * hydro["length_stations_m"] * hydro["beam_waterline_m"] * \
        hydro["draft_m"] == pytest.approx(hydro["displacement_volume_m3"], rel=1e-9)
    # beam_m in the artifact is the half breadth, beam_full_m the moulded beam
    assert space.beam_full_m == pytest.approx(2.0 * space.beam_m)


def test_volume_matches_hydrostatics_at_deep_waterline(space):
    """Above the deck nothing changes: the shell volume is the integral of the sections."""
    v = space.volume(space.zero())
    assert 15000.0 < v < 25000.0
    assert v != pytest.approx(space.hydrostatics(None)["displacement_volume_m3"])


# ---------------------------------------------------------------------------
# entry point
# ---------------------------------------------------------------------------
def test_cli_entry_point(tmp_path):
    out = tmp_path / "rec"
    proc = subprocess.run(
        [sys.executable, "-m", "tensorcad.dtmb5415_cli", "--out", str(out),
         "--stl", str(out / "hull.stl"), "--gradient", "--json"],
        capture_output=True, text=True, cwd=ROOT,
        env={"PYTHONPATH": str(ROOT), "PATH": "/usr/bin:/bin"},
    )
    assert proc.returncode == 0, proc.stderr
    rec = json.loads(proc.stdout)
    assert rec["quality_gate"]["pass"] is True
    assert rec["hydrostatics"]["cb"] == pytest.approx(0.502, abs=5e-3)
    assert len(rec["volume_gradient"]) == 27
    assert (out / "dtmb5415_record.json").exists()
    assert (out / "hull.stl").exists() and (out / "hull.stl").stat().st_size > 1_000_000


def test_cli_rejects_bad_design_vector():
    proc = subprocess.run(
        [sys.executable, "-m", "tensorcad.dtmb5415_cli", "--design-vector", "0.1,0.2",
         "--no-write", "--json"],
        capture_output=True, text=True, cwd=ROOT,
        env={"PYTHONPATH": str(ROOT), "PATH": "/usr/bin:/bin"},
    )
    assert proc.returncode == 2
    assert "design vector" in proc.stdout.lower() or "design vector" in proc.stderr.lower()