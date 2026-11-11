"""Tests for the parametric hull command-line entry point.

The entry point is the only way a user (or an agent) reaches the hull pipeline without
writing Python, so it is tested the way it is used: through ``main()`` with argument
lists, checking the JSON record it emits and the gmsh file it writes.
"""
from __future__ import annotations

import json

import pytest

from tensorcad.adapters.fem import read_gmsh
from tensorcad.parametric_hull_cli import main


def _run(argv, capsys):
    code = main(argv)
    captured = capsys.readouterr()
    return code, captured.out, captured.err


def test_wigley_run_reports_analytic_comparison(tmp_path, capsys):
    record = tmp_path / "wigley.json"
    code, out, _ = _run(
        ["wigley", "--length", "10", "--beam", "1", "--draft", "0.625",
         "--json", str(record)],
        capsys,
    )
    assert code == 0
    payload = json.loads(out)
    assert json.loads(record.read_text()) == payload
    assert payload["gate_failures"] == []
    comparison = payload["analytic_comparison"]
    # Criterion (b) of the node: displacement within 0.5 % of the closed form.
    assert abs(comparison["volume_rel_error"]) < 5e-3
    # Criterion (c): waterplane area within 1 %.
    assert abs(comparison["waterplane_area_rel_error"]) < 1e-2
    assert payload["hydrostatics"]["displaced_volume_m3"] > 0.0


def test_wigley_cli_honours_requested_block_coefficient(tmp_path, capsys):
    code, out, _ = _run(
        ["wigley", "--block-coefficient", "0.55", "--mesh-n-u", "41", "--mesh-n-v", "10",
         "--json", str(tmp_path / "cb.json")],
        capsys,
    )
    assert code == 0
    payload = json.loads(out)
    assert payload["definition"]["requested_block_coefficient"] == pytest.approx(0.55)
    assert payload["hydrostatics"]["block_coefficient"] == pytest.approx(0.55, abs=5e-3)


def test_cli_gmsh_handoff_is_watertight_and_closed(tmp_path, capsys):
    msh = tmp_path / "hull.msh"
    code, out, _ = _run(
        ["wigley", "--length", "10", "--beam", "1", "--draft", "0.625",
         "--mesh-n-u", "41", "--mesh-n-v", "10", "--gmsh", str(msh)],
        capsys,
    )
    assert code == 0
    payload = json.loads(out)
    assert payload["gmsh"]["watertight"] is True
    assert payload["gmsh"]["edge_manifold"] is True
    # A closed hull is a topological sphere: V - E + F = 2.
    assert payload["gmsh"]["euler_characteristic"] == 2
    parsed = read_gmsh(msh)
    assert parsed["nodes"] == payload["mesh"]["vertices"]
    assert parsed["triangles"] == payload["mesh"]["faces"]
    # The volume the file carries is the volume the hydrostatics report (same mesh).
    assert payload["gmsh"]["volume_m3"] == pytest.approx(
        payload["hydrostatics"]["displaced_volume_m3"], rel=1e-12
    )


def test_cli_suboff_with_clustered_stations(tmp_path, capsys):
    code, out, _ = _run(
        ["suboff", "--station-fractions", "--stations", "21",
         "--json", str(tmp_path / "suboff.json")],
        capsys,
    )
    assert code == 0
    payload = json.loads(out)
    assert payload["definition"]["clustered_stations"] is True
    assert payload["sections"] > 21  # the pinned junctions and stern-cap stations are added
    assert payload["hydrostatics"]["displaced_volume_m3"] > 0.0
    assert payload["gate_failures"] == []
    # A body of revolution has no coplanar section tops, so the lid is skipped and the
    # record states it rather than silently producing a wrong closure.
    assert "waterplane_lid" in payload["definition"]
    assert payload["mesh"]["waterplane_closed"] is False


def test_cli_suboff_gmsh_handoff_is_refused_not_wrong(tmp_path, capsys):
    """An open shell cannot be exported as a body; the CLI must say so, not guess."""
    code, _out, err = _run(
        ["suboff", "--stations", "11", "--mesh-n-u", "21", "--mesh-n-v", "8",
         "--gmsh", str(tmp_path / "suboff.msh")],
        capsys,
    )
    assert code == 2
    assert "waterplane lid" in err


def test_cli_offset_table_reproduces_wedge_volume(tmp_path, capsys):
    """A wedge section (y grows linearly from the keel) lofts into V = L * b_wl * T exactly."""
    length, half_beam, draft = 10.0, 0.5, 0.5
    rows = 8
    points = [
        [half_beam * i / (rows - 1), -draft + draft * i / (rows - 1)] for i in range(rows)
    ]
    table = {
        "length": length,
        "stations": [
            {"x_fraction": frac, "points": points} for frac in (0.0, 0.25, 0.5, 0.75, 1.0)
        ],
    }
    path = tmp_path / "wedge.json"
    path.write_text(json.dumps(table))
    code, out, _ = _run(
        ["offsets", str(path), "--mesh-n-u", "41", "--mesh-n-v", "6",
         "--json", str(tmp_path / "wedge_out.json")],
        capsys,
    )
    assert code == 0
    payload = json.loads(out)
    volume = payload["hydrostatics"]["displaced_volume_m3"]
    assert volume == pytest.approx(length * half_beam * draft, rel=1e-6)


def test_cli_rejects_open_waterplane_handoff(tmp_path, capsys):
    code, _, err = _run(
        ["wigley", "--mesh-n-u", "21", "--mesh-n-v", "8", "--open-waterplane",
         "--gmsh", str(tmp_path / "open.msh")],
        capsys,
    )
    assert code == 2
    assert "close_waterplane" in err


def test_cli_missing_offset_table_is_a_usage_error(tmp_path, capsys):
    code, _, err = _run(["offsets", str(tmp_path / "nope.json")], capsys)
    assert code == 2
    assert "does not exist" in err


def test_parser_requires_a_family():
    with pytest.raises(SystemExit):
        main([])


def test_cli_offset_table_accepts_absolute_and_fractional_stations(tmp_path, capsys):
    def block(x_key, x_value):
        return {
            x_key: x_value,
            "half_breadths": [0.0, 0.25, 0.5],
            "heights": [-0.5, -0.25, 0.0],
        }

    table = {"length": 8.0, "stations": [block("x_fraction", 0.0), block("x", 8.0)]}
    path = tmp_path / "mixed.json"
    path.write_text(json.dumps(table))
    code, out, _ = _run(
        ["offsets", str(path), "--mesh-n-u", "21", "--mesh-n-v", "6",
         "--json", str(tmp_path / "mixed_out.json")],
        capsys,
    )
    assert code == 0
    payload = json.loads(out)
    assert payload["definition"]["family"] == "offset_table"
    assert payload["sections"] == 2
    # A V-section over L = 8 m with half-breadth 0.5 m at the waterline and 0.5 m draft:
    # V = L * b_wl * T (both sides) = 2.0 m^3.
    assert payload["hydrostatics"]["displaced_volume_m3"] == pytest.approx(2.0, rel=2e-3)


@pytest.mark.parametrize(
    ("payload", "message"),
    [
        ({"stations": []}, "no 'stations' list"),
        ({"stations": [{"points": [[0.1, -0.5], [0.2, 0.0]]}]}, "neither 'x' nor 'x_fraction'"),
        ({"stations": [{"x": 0.0, "points": [[0.1, -0.5]]}]}, "at least 2 points"),
    ],
)
def test_cli_offset_table_rejects_malformed_input(tmp_path, capsys, payload, message):
    path = tmp_path / "bad.json"
    path.write_text(json.dumps(payload))
    code, _, err = _run(["offsets", str(path)], capsys)
    assert code == 2
    assert message in err


def test_cli_gate_failures_are_reported_not_swallowed(monkeypatch, capsys):
    """A degenerate build must still emit its record, with the failure named and rc=1."""
    from tensorcad import parametric_hull_cli as cli

    real_build = cli.build_parametric_hull

    def broken(*args, **kwargs):
        hull = real_build(*args, **kwargs)
        hull.quality = dict(hull.quality, min_scaled_jacobian=0.0)
        return hull

    monkeypatch.setattr(cli, "build_parametric_hull", broken)
    code, out, _ = _run(["wigley", "--mesh-n-u", "21", "--mesh-n-v", "8"], capsys)
    assert code == 1
    payload = json.loads(out)
    assert payload["gate_failures"]
    assert "min_scaled_jacobian" in payload["gate_failures"][0]
