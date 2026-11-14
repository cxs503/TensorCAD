"""Tests for adopting the derived JBC model scale (anchor GAP-CAD-T1-SOURCE-HUNT ③).

The derivation itself is tested in ``test_jbc_model_scale.py``.  What is pinned here is
the *plumbing* that was missing: the record must be committed (a skipped record test is
not evidence), the fleet manifest must use it, and the production route must not fall
back to the ``L_PP`` x-extent fit that shrinks the hull by ~5 %.

Nothing here rebuilds the hull: the committed impact record carries the two build
reports, and these tests assert its invariants offline.
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
RECORD = ROOT / "artifacts" / "baseline_fleet" / "jbc_model_scale.json"
IMPACT = ROOT / "artifacts" / "baseline_fleet" / "jbc_scale_impact.json"
SPECS = ROOT / "benchmarks" / "verified" / "baseline_fleet" / "fleet_specs.json"
ADOPT = ROOT / "tools" / "baseline_fleet" / "adopt_derived_scale.py"
BUILD = ROOT / "tools" / "baseline_fleet" / "build_fleet.py"


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.path.insert(0, str(ROOT))
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def record():
    assert RECORD.exists(), (
        "the derived-scale record is committed evidence: without it the record tests in "
        "test_jbc_model_scale.py skip and the published numbers are unreproducible"
    )
    return json.loads(RECORD.read_text())


def test_record_is_schema_tagged_and_digest_bound(record):
    assert record["schema"] == "tensorcad.fleet.model_scale.v1"
    assert record["anchor"] == "GAP-CAD-T1-SOURCE-HUNT"
    assert record["scale_derived"]["value"] == pytest.approx(39.72849598755336, rel=1e-12)
    assert record["scale_derived"]["relative_half_spread"] < 0.02
    for key in ("particulars_page", "iges_archive"):
        assert len(record["sources"][key]["sha256"]) == 64
    assert record["checks"]["iges_header_units"]["unit_name"] == "MM"


def test_manifest_uses_the_derived_scale_not_the_lpp_fit(record):
    ships = {ship["name"]: ship for ship in json.loads(SPECS.read_text())["ships"]}
    jbc = ships["JBC"]
    assert jbc["geometry_scale"] is not None, (
        "geometry_scale None means the builder fits the x-extent to L_PP, which builds "
        "JBC ~5 % short"
    )
    assert jbc["geometry_scale"] == pytest.approx(record["scale_derived"]["value"], rel=1e-12)
    assert "jbc_model_scale.json" in jbc["scale_basis"]
    # the other three hulls must keep their own basis
    assert ships["KCS"]["geometry_scale"] == 1.0
    assert ships["KVLCC2"]["geometry_scale"] == 1.0
    assert ships["ONRT"]["geometry_scale"] == pytest.approx(48.9355, rel=1e-6)


def test_adopt_tool_check_passes_on_the_committed_manifest():
    module = _load(ADOPT, "adopt_derived_scale")
    assert module.sync(SPECS, RECORD, write=False) == 0


def test_adopt_tool_is_idempotent_and_refuses_unknown_schema(tmp_path, record):
    module = _load(ADOPT, "adopt_derived_scale_idem")
    specs = tmp_path / "fleet_specs.json"
    specs.write_text(json.dumps({"ships": [{"name": "JBC", "geometry_scale": None}]}))
    assert module.sync(specs, RECORD, write=True) == 0
    first = specs.read_text()
    assert module.sync(specs, RECORD, write=True) == 0
    assert specs.read_text() == first, "second run must not churn the manifest"

    bad = tmp_path / "bad.json"
    bad.write_text(json.dumps({"schema": "something.else", "scale_derived": {"value": 1.0}}))
    with pytest.raises(SystemExit):
        module.derived_scale(bad)


def test_build_policy_routes_jbc_through_the_record():
    module = _load(BUILD, "build_fleet")
    assert module.GEOMETRY_SCALE_POLICY["JBC"] == ("derived", "jbc")
    derived = module.derived_scale("jbc")
    assert derived is not None
    scale, basis = derived
    assert scale == pytest.approx(39.72849598755336, rel=1e-12)
    assert "jbc_model_scale.json" in basis and "sha256" in basis
    assert module.derived_scale("no-such-hull") is None


def test_build_policy_falls_back_without_the_record():
    module = _load(BUILD, "build_fleet_fallback")
    module.DERIVED_SCALE_RECORDS["jbc"] = ROOT / "artifacts" / "baseline_fleet" / "missing.json"
    try:
        assert module.derived_scale("jbc") is None
    finally:
        module.DERIVED_SCALE_RECORDS["jbc"] = RECORD


def test_impact_record_quantifies_the_bias_correctly():
    data = json.loads(IMPACT.read_text())
    inferred, derived = data["variants"]["inferred"], data["variants"]["derived"]
    # the x-extent fit is the biased reference: it shrinks the hull
    assert inferred["scale_inferred"] is True and derived["scale_inferred"] is False
    assert inferred["length_waterline_m"] < derived["length_waterline_m"]
    assert -10.0 < inferred["error_percent"]["length_waterline_m"] < -3.0
    assert -10.0 < inferred["error_percent"]["beam_waterline_m"] < -3.0
    # adopting the scale brings the lengths inside the 3 % acceptance band
    assert abs(derived["error_percent"]["length_waterline_m"]) < 3.0
    assert abs(derived["error_percent"]["beam_waterline_m"]) < 3.0
    # the displacement deficit does NOT go away: the publisher archive is not a closed
    # shell, so this stays an honest FAIL and must not be presented as fixed
    assert derived["error_percent"]["volume_m3"] < -5.0
    assert derived["verdict"] == "WITHHELD"