"""Tests for the real-CAD-source performance ladder (GAP-CAD-D4-PERF-SCALE).

The load-bearing assertions are:

* the **ULP domain defect** that made a real IGES hull unreadable stays fixed
  (``test_patch_boundary_sampler_stays_inside_the_domain`` +
  ``test_onrt_patch_set_reads_after_the_endpoint_fix``);
* **chunked tessellation is bit-identical** to the one-shot mesh, and the gate
  that says so can actually fail (a gate whose RED path is never exercised is
  not a gate);
* the **2x wall-clock gate** turns RED at 3x and reports ``not_run`` - not a
  pass - when the baseline has no key;
* a rung above ``max_triangles`` is ``not_run`` with a reason instead of being
  silently skipped or hard-run into an OOM.
"""
from __future__ import annotations

import math
from pathlib import Path

import pytest

from tensorcad import perf_cad_source as pcs
from tensorcad.geomquality.brep import _domain_param

REPO = Path(__file__).resolve().parents[1]
ONRT = REPO / "work" / "assets" / "onrt.igs"


# ---------------------------------------------------------------------------
# the ULP domain defect
# ---------------------------------------------------------------------------
def test_domain_param_pins_the_endpoints_exactly():
    """``a + (b - a) * t`` can leave the domain at t == 1; the helper may not."""
    a, b = -35.07178637012983, 2.619396642977591
    # Precondition: the naive interpolation really does overshoot this domain,
    # so the test is testing a live hazard rather than a hypothetical one.
    naive = a + (b - a) * 1.0
    assert naive > b
    assert _domain_param(a, b, 1.0) == b
    assert _domain_param(a, b, 0.0) == a
    # Interior values are untouched (no silent clamping of legitimate samples).
    mid = _domain_param(a, b, 0.5)
    assert a < mid < b
    assert mid == a + (b - a) * 0.5


def test_domain_param_clamps_an_out_of_range_fraction():
    a, b = 0.0, 3.0
    assert _domain_param(a, b, 1.5) == b
    assert _domain_param(a, b, -0.5) == a


def test_onrt_patch_set_reads_after_the_endpoint_fix():
    """A tracked real IGES hull must be readable, with a usable topology."""
    if not ONRT.exists():  # pragma: no cover - source missing in a slim checkout
        pytest.skip("tracked ONRT sample not present")
    brep, seconds, census = pcs.read_source(ONRT, pcs.ROUTE_PATCH)
    assert census["faces"] == 132
    assert census["edges"] > 0
    assert seconds >= 0.0
    assert brep.faces


# ---------------------------------------------------------------------------
# chunked tessellation
# ---------------------------------------------------------------------------
@pytest.fixture(scope="module")
def onrt_brep():
    if not ONRT.exists():  # pragma: no cover
        pytest.skip("tracked ONRT sample not present")
    return pcs.read_source(ONRT, pcs.ROUTE_PATCH)[0]


def test_chunked_tessellation_is_bit_identical(onrt_brep):
    reference = onrt_brep.tessellate(4, 3)
    assert len(reference.faces) > 0
    for chunk in (1, 5, 64, None):
        chunked = onrt_brep.tessellate(4, 3, chunk_faces=chunk)
        assert tuple(chunked.vertices) == tuple(reference.vertices)
        assert tuple(chunked.faces) == tuple(reference.faces)


def test_chunk_callback_visits_every_face_once(onrt_brep):
    seen: list[tuple[int, int]] = []
    onrt_brep.tessellate(2, 2, chunk_faces=8, on_chunk=lambda done, total: seen.append((done, total)))
    assert seen
    assert seen[-1][0] == seen[-1][1] == len(onrt_brep.faces)
    assert [done for done, _ in seen] == sorted({done for done, _ in seen})


def test_chunk_faces_must_be_positive(onrt_brep):
    with pytest.raises(ValueError):
        onrt_brep.tessellate(2, 2, chunk_faces=0)


def test_chunk_equivalence_gate_fires_red_on_a_changed_mesh():
    good = {"route": "patch", "parameter": "uv=4x3", "chunked": {"identical": True}}
    bad = {"route": "patch", "parameter": "uv=4x3", "chunked": {"identical": False}}
    gates = pcs.chunk_equivalence_gates([good, bad])
    assert [g.passed for g in gates] == [True, False]
    assert "RED" in gates[1].detail
    # A rung that was never chunked contributes no verdict at all.
    assert pcs.chunk_equivalence_gates([{"route": "patch", "parameter": "uv=4x3"}]) == ()


# ---------------------------------------------------------------------------
# the 2x regression gate
# ---------------------------------------------------------------------------
def _record(seconds=0.1, weld=0.02, volume=0.01, locate=0.05, chunked=0.11):
    return {
        "route": "patch",
        "parameter": "uv=8x6",
        "status": pcs.MEASURED,
        "seconds": seconds,
        "weld": {"seconds": weld},
        "volume": {"seconds": volume},
        "locate": {"seconds": locate},
        "chunked": {"seconds": chunked, "identical": True},
    }


def test_regression_gate_is_red_at_three_times_the_baseline():
    baseline = {"measurements": {
        pcs.measurement_key("cad_tessellate", "deadbeef" * 8, "patch", "uv=8x6"): 0.1,
        pcs.measurement_key("cad_weld", "deadbeef" * 8, "patch", "uv=8x6"): 0.02,
        pcs.measurement_key("cad_volume", "deadbeef" * 8, "patch", "uv=8x6"): 0.01,
        pcs.measurement_key("cad_locate", "deadbeef" * 8, "patch", "uv=8x6"): 0.05,
        pcs.measurement_key("cad_tessellate_chunked", "deadbeef" * 8, "patch", "uv=8x6"): 0.11,
    }}
    gates = pcs.regression_gates([_record(seconds=0.3)], baseline, "deadbeef" * 8)
    tessellation = [g for g in gates if g.name.startswith("regression_cad_tessellate[")][0]
    assert tessellation.passed is False
    assert "RED" in tessellation.detail
    # The other sub-tasks were unchanged and stay green.
    assert all(g.passed for g in gates if g is not tessellation)


def test_regression_gate_passes_below_the_threshold():
    key = pcs.measurement_key("cad_tessellate", "cafe" * 16, "patch", "uv=8x6")
    gates = pcs.regression_gates([_record(seconds=0.19)], {"measurements": {key: 0.1}}, "cafe" * 16)
    tessellation = [g for g in gates if g.name.startswith("regression_cad_tessellate[")][0]
    assert tessellation.passed is True


def test_regression_gate_without_a_baseline_key_is_not_run_not_pass():
    gates = pcs.regression_gates([_record()], {"measurements": {}}, "beef" * 16)
    assert gates
    assert all(g.status == pcs.NOT_RUN for g in gates)
    # ``passed`` is True for not_run (backwards compatible) but the *status*
    # says "not checked"; the ladder verdict must not read that as a pass.
    assert pcs._verdict(gates) == "not_run"


def test_regression_gate_ignores_a_rung_that_was_not_measured():
    record = {"route": "patch", "parameter": "uv=8x6", "status": pcs.NOT_RUN,
              "seconds": None, "note": "not_run: too big"}
    gates = pcs.regression_gates([record], {"measurements": {}}, "abcd" * 16)
    assert all(g.status == pcs.NOT_RUN for g in gates)
    assert any("no measurement" in g.detail for g in gates)


# ---------------------------------------------------------------------------
# scale safety and cluster thresholds
# ---------------------------------------------------------------------------
def test_a_rung_above_max_triangles_is_not_run_with_a_reason(onrt_brep, tmp_path):
    rung = pcs.Rung(u_samples=400, v_samples=400)  # 132 * 2 * 400 * 400 = 42.2e6
    record, mesh = pcs.measure_rung(ONRT, pcs.ROUTE_PATCH, rung, model=onrt_brep,
                                    max_triangles=16_000_000)
    assert mesh is None
    assert record["status"] == pcs.NOT_RUN
    assert record["seconds"] is None
    assert "max_triangles" in record["note"]
    assert "42,240,000" in record["note"]


def test_cluster_threshold_needs_two_points():
    out = pcs.cluster_thresholds([{"status": pcs.MEASURED, "triangles": 1000, "seconds": 0.1,
                                   "route": "patch", "parameter": "uv=4x3"}])
    assert out["cad_tessellate"]["status"] == pcs.NOT_RUN
    assert "2 measured sizes" in out["cad_tessellate"]["reason"]


def test_cluster_threshold_inverts_a_linear_fit():
    records = [
        {"status": pcs.MEASURED, "triangles": 1_000, "seconds": 0.001,
         "route": "patch", "parameter": "uv=1x1"},
        {"status": pcs.MEASURED, "triangles": 1_000_000, "seconds": 1.0,
         "route": "patch", "parameter": "uv=2x2"},
    ]
    out = pcs.cluster_thresholds(records, budget_seconds=60.0)["cad_tessellate"]
    assert out["exponent"] == pytest.approx(1.0, abs=1e-6)
    # 60x the 1 s / 1e6-triangle rate, exponent 1 -> 6.0e7 triangles.
    assert out["size_at_budget"] == pytest.approx(6.0e7, rel=1e-3)
    assert out["exceeds_max_rung"] is True


# ---------------------------------------------------------------------------
# rung plumbing
# ---------------------------------------------------------------------------
def test_parse_rungs_handles_both_axes():
    patch = pcs.parse_rungs("8x6, 16x12")
    assert [(r.u_samples, r.v_samples) for r in patch] == [(8, 6), (16, 12)]
    trimmed = pcs.parse_rungs("chord=0.01,chord=0.001")
    assert [r.chord_tolerance for r in trimmed] == [0.01, 0.001]
    with pytest.raises(ValueError):
        pcs.parse_rungs("   ")


def test_parameter_and_expected_triangles():
    rung = pcs.Rung(u_samples=4, v_samples=3)
    assert rung.parameter(pcs.ROUTE_PATCH) == "uv=4x3"
    assert rung.expected_triangles(132) == 2 * 4 * 3 * 132
    trimmed = pcs.Rung(chord_tolerance=0.01)
    assert trimmed.parameter(pcs.ROUTE_TRIMMED) == "chord=0.01"
    with pytest.raises(ValueError):
        pcs.Rung(u_samples=1, v_samples=1).parameter(pcs.ROUTE_TRIMMED)


def test_measurement_key_separates_parameters_and_sources():
    a = pcs.measurement_key("cad_weld", "aa" * 32, "patch", "uv=8x6")
    b = pcs.measurement_key("cad_weld", "aa" * 32, "patch", "uv=16x12")
    c = pcs.measurement_key("cad_weld", "bb" * 32, "patch", "uv=8x6")
    assert len({a, b, c}) == 3
    assert a.split("|")[1] == "a" * 16  # digest prefix is in the key


def test_full_ladder_smoke_on_the_real_source(tmp_path):
    """End to end on the tracked source: measure, gate, save, re-gate."""
    if not ONRT.exists():  # pragma: no cover
        pytest.skip("tracked ONRT sample not present")
    payload = pcs.run_cad_ladder(ONRT, (pcs.Rung(4, 3), pcs.Rung(8, 6)),
                                 route=pcs.ROUTE_PATCH, queries=16, chunk_faces=8,
                                 max_triangles=1_000_000)
    assert payload["schema"] == pcs.SCHEMA
    assert len(payload["records"]) == 2
    first, second = payload["records"]
    assert first["triangles"] > 0
    assert second["triangles"] > first["triangles"]
    assert first["chunked"]["identical"] is True
    weld = first["weld"]
    assert weld["vertices_merged"] >= 0
    assert weld["boundary_edges_after"] <= weld["boundary_edges_before"]
    assert math.isfinite(first["volume"]["signed_volume"])
    assert first["locate"]["queries"] == 16
    # No baseline yet -> the regression gates are not_run, so the verdict is
    # not_run (explicitly "not checked") rather than a green pass.
    assert payload["verdict"] == "not_run"
    assert payload["gates_not_run"] > 0

    # Round trip the baseline and re-gate against it: now the gates are live.
    target = tmp_path / "baseline.json"
    target.write_text(__import__("json").dumps(
        {"schema": pcs.SCHEMA, "measurements": payload["measurements"]}))
    import json

    again = pcs.run_cad_ladder(ONRT, (pcs.Rung(4, 3),), route=pcs.ROUTE_PATCH,
                               queries=16, chunk_faces=8, max_triangles=1_000_000,
                               baseline=json.loads(target.read_text()))
    regression = [g for g in again["gates"] if g["name"].startswith("regression_cad_tessellate[")]
    assert regression and regression[0]["status"] == pcs.MEASURED
    assert again["verdict"] in ("pass", "fail")