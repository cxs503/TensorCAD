"""Tests for the JBC model-scale derivation (anchor GAP-CAD-T1-SOURCE-HUNT item 3).

Two flavours of check:

* **arithmetic / parser** tests run on synthetic inputs, so they are hermetic and pin
  the invariants the derivation relies on (span -> scale, the two-row table layout,
  the ``(m 3 )`` unit trap that once turned the displacement into ``3.0``);
* **record** tests re-derive the committed numbers offline.  They never download:
  archive bytes are deliberately not redistributed, only digests and spans are.
"""

from __future__ import annotations

import importlib.util
import json
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "tools" / "baseline_fleet" / "jbc_model_scale.py"
RECORD_PATH = ROOT / "artifacts" / "baseline_fleet" / "jbc_model_scale.json"


def _load():
    spec = importlib.util.spec_from_file_location("jbc_model_scale", MODULE_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def mod():
    return _load()


class StubSurface:
    """A planar patch: evaluate(u, v) -> (a*u, b*v, c). Enough to pin surface_bbox."""

    def __init__(self, a: float = 2.0, b: float = 3.0, c: float = 0.0):
        self.a, self.b, self.c = a, b, c
        self.domain_u = (0.0, 1.0)
        self.domain_v = (0.0, 1.0)

    def evaluate(self, u: float, v: float):
        return (self.a * u, self.b * v, self.c)


def test_surface_bbox_measures_the_sampled_extent(mod):
    box = mod.surface_bbox([StubSurface(), StubSurface(c=5.0)], samples=11)
    assert box["surfaces"] == 2
    assert box["points"] == 2 * 11 * 11
    assert box["span"] == pytest.approx([2.0, 3.0, 5.0])
    assert box["lo"][2] == pytest.approx(0.0)
    assert box["hi"][2] == pytest.approx(5.0)


def test_surface_bbox_survives_a_degenerate_corner(mod):
    class Flaky(StubSurface):
        def evaluate(self, u, v):
            if (u, v) == (0.0, 0.0):
                raise ValueError("singular corner")
            return super().evaluate(u, v)

    box = mod.surface_bbox([Flaky()], samples=5)
    assert box["points"] == 5 * 5 - 1
    assert box["span"] == pytest.approx([2.0, 3.0, 0.0])


def test_derive_scale_recovers_a_known_scale(mod):
    published = {"beam_m": 45.0, "depth_m": 25.0, "lpp_m": 280.0}
    # A hull drawn at exactly 1:40
    span = {"x": 7.5, "y": 45.0 / 40.0, "z": 25.0 / 40.0}
    scale = mod.derive_scale(span, published)
    assert scale["value"] == pytest.approx(40.0)
    assert scale["relative_half_spread"] == pytest.approx(0.0, abs=1e-12)
    assert abs(scale["hypothesis_residuals"]["1:40"]["beam"]) < 1e-12
    # the x extent is recorded but must not be used as an Lpp reference
    assert scale["not_used"]["x_span"]["naive_scale_against_lpp"] == pytest.approx(280.0 / 7.5)


def test_derive_scale_reports_the_reference_disagreement(mod):
    published = {"beam_m": 45.0, "depth_m": 25.0, "lpp_m": 280.0}
    # depth 1 % short of the beam-derived scale: that disagreement must show up
    span = {"x": 7.5, "y": 45.0 / 39.8, "z": 25.0 / 39.4}
    scale = mod.derive_scale(span, published)
    assert scale["relative_half_spread"] > 0.004
    assert scale["candidates"]["beam"]["scale"] == pytest.approx(39.8)
    assert scale["candidates"]["depth"]["scale"] == pytest.approx(39.4)


def test_particulars_parser_ignores_unit_digits(mod):
    """Regression: the label row carries digits inside its unit ("(m 3 )")."""
    html = """
    <table>
      <tr><td>Length between perpendiculars</td><td>L<sub>PP</sub> (m)</td></tr>
      <tr><td>280.000</td><td>same as on the left</td></tr>
      <tr><td>Maximum beam of waterline</td><td>B<sub>WL</sub> (m)</td></tr>
      <tr><td>45.000</td><td>same as on the left</td></tr>
      <tr><td>Depth</td><td>D (m)</td></tr>
      <tr><td>25.000</td><td>same as on the left</td></tr>
      <tr><td>Displacement volume of a hull</td><td>&nabla; (m<sup>3</sup>)</td></tr>
      <tr><td>178369.9</td><td>89185.9</td></tr>
    </table>
    """
    out = mod.read_particulars(html)
    assert out["lpp_m"]["value"] == 280.0
    assert out["beam_m"]["value"] == 45.0
    assert out["depth_m"]["value"] == 25.0
    assert out["displacement_m3"]["value"] == 178369.9   # not 3.0 from "(m 3 )"
    assert "178369.9" in out["displacement_m3"]["snippet"]


def test_model_propeller_parser(mod):
    text = "PROPELLER NO.                MP. 687\nDIAMETER                      203.00 MM\n"
    out = mod.read_model_propeller(text)
    assert out["diameter_mm"] == 203.0
    assert "DIAMETER" in out["snippet"]


@pytest.mark.skipif(not RECORD_PATH.exists(), reason="record not generated")
def test_committed_record_verifies_offline(mod, capsys):
    assert mod.verify(RECORD_PATH) == 0
    assert "OK" in capsys.readouterr().out


@pytest.mark.skipif(not RECORD_PATH.exists(), reason="record not generated")
def test_committed_record_binds_every_number_to_a_source():
    record = json.loads(RECORD_PATH.read_text())
    assert record["anchor"] == "GAP-CAD-T1-SOURCE-HUNT"
    for name in ("particulars_page", "iges_archive"):
        info = record["sources"][name]
        assert re.fullmatch(r"[0-9a-f]{64}", info["sha256"]), name
        assert info["http_status"] == 200 and info["bytes"] > 0
    # every published value arrives with the row text it was read from
    for key, item in record["published_full_scale"].items():
        assert item["snippet"], key
        assert item["value"] is not None, key
    assert record["published_full_scale"]["displacement_m3"]["value"] == pytest.approx(178369.9)
    assert record["published_full_scale"]["lpp_m"]["value"] == pytest.approx(280.0)

    geom = record["model_geometry"]
    assert geom["surfaces"] > 10 and geom["sampled_points"] > 10000
    assert record["checks"]["iges_header_units"]["unit_name"] == "MM"
    assert record["checks"]["keel_plane_at_zero"]["within_1mm"] is True
    assert record["checks"]["half_model_symmetric"]["asymmetry_m"] < 1e-9

    scale = record["scale_derived"]
    # two independent references, agreeing within 2 % (the verify() bound)
    assert len(scale["candidates"]) == 2
    assert 0.0 < scale["relative_half_spread"] < 0.02
    # the derived scale must reproduce the published Cb through its own model volume
    implied = record["checks"]["implied_model_lengths"]
    cb = record["published_full_scale"]["cb"]["value"]
    closure = (record["published_full_scale"]["displacement_m3"]["value"]
               / (record["published_full_scale"]["lpp_m"]["value"]
                  * record["published_full_scale"]["beam_m"]["value"]
                  * record["published_full_scale"]["draft_m"]["value"]))
    assert closure == pytest.approx(cb, rel=1e-3)
    assert implied["lpp_m"] == pytest.approx(
        record["published_full_scale"]["lpp_m"]["value"] / scale["value"])
    # the hull's x extent must exceed the implied Lpp (deck-level overhang), which is
    # exactly why it is excluded as a scale reference
    assert implied["x_extent_over_implied_lpp"] > 1.0