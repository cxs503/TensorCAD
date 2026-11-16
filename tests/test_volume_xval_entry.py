"""Entry and integration tests for the cross-tool volume cross-check.

Three things are checked here that the engine tests in ``test_volume_xval.py``
cannot check, because they are properties of how the capability is *reached*:

* the package surface exports it (a caller who does not know the module path can
  still find it), and the console script dispatches to it;
* an operator gets a *refusal* (exit 2) and an *honest FAIL* (exit 1) as well as a
  pass (exit 0), with the JSON evidence written either way;
* the product gate that grades hull volumes -- :func:`tensorcad.fleet_gate.audit_fleet`,
  the gate the JBC ``-18.3 %`` episode went through -- actually calls the independent
  engines, records what it measured, and refuses to certify a volume the engines
  cannot reproduce.
"""

from __future__ import annotations

import json
import struct
from pathlib import Path

import pytest

import tensorcad
from tensorcad import volume_xval as xv
from tensorcad.fleet_gate import VERIFIED, audit_fleet
from tensorcad.geometry import Mesh

EXPORTED = (
    "ACCEPTANCE_TOLERANCE",
    "ENGINE_AGREEMENT_TOLERANCE",
    "VolumeXvalError",
    "classify_routes",
    "closed_solid_from_breadth",
    "engine_agreement_note",
    "fit_order",
    "richardson_limit",
    "topology_report",
    "volume_divergence",
    "volume_numpy_stl",
    "volume_sectional",
    "volume_trimesh",
)

trimesh = pytest.importorskip("trimesh", reason="the independent engine route needs trimesh")


def _write_binary_stl(path: Path, mesh: Mesh) -> Path:
    with path.open("wb") as handle:
        handle.write(b"tensorcad volume_xval entry test".ljust(80, b" "))
        handle.write(struct.pack("<I", len(mesh.faces)))
        for face in mesh.faces:
            handle.write(struct.pack("<3f", 0.0, 0.0, 0.0))
            for index in face:
                handle.write(struct.pack("<3f", *mesh.vertices[index]))
            handle.write(struct.pack("<H", 0))  # per-facet attribute byte count
    return path


def _box_stl(path: Path, extents=(2.0, 3.0, 5.0)) -> Path:
    box = trimesh.creation.box(extents=extents)
    mesh = Mesh(vertices=[tuple(map(float, v)) for v in box.vertices], faces=[tuple(map(int, f)) for f in box.faces])
    return _write_binary_stl(path, mesh)


def _wigley_solid(*, stations: int = 81, levels: int = 17) -> tuple[Mesh, float]:
    from tensorcad import wigley as wg

    parameters = wg.reference_parameters()
    solid = xv.closed_solid_from_breadth(
        lambda x, z: wg.half_breadth(parameters, x, z),
        length=parameters.length,
        draft=parameters.draft,
        stations=stations,
        levels=levels,
    )
    return solid, wg.analytic_volume(parameters.length, parameters.beam, parameters.draft)


# ---------------------------------------------------------------------------
# package surface
# ---------------------------------------------------------------------------
def test_package_exports_the_cross_check():
    for name in EXPORTED:
        assert name in tensorcad.__all__, name
        assert hasattr(tensorcad, name), name
    assert tensorcad.ACCEPTANCE_TOLERANCE == xv.ACCEPTANCE_TOLERANCE


def test_classifier_refuses_to_name_a_bias_without_a_ladder():
    """A gap with no measured order is unattributed -- not 'systematic bias'."""
    verdict = xv.classify_routes(divergence=1.0, engine_a=1.0, engine_b=1.0, pipeline=1.01)
    assert verdict["status"] == "undetermined"
    assert "no refinement ladder" in verdict["reason"]
    assert verdict["observed_order"] is None


# ---------------------------------------------------------------------------
# console entry
# ---------------------------------------------------------------------------
def test_console_script_reaches_the_cross_check(capsys):
    from tensorcad.cli import main as cli_main

    code = cli_main(
        ["volume-xval", "wigley", "--stations", "81", "--levels", "17", "--json"]
    )
    payload = json.loads(capsys.readouterr().out)
    assert code == 0
    assert payload["case"]["case"] == "wigley"
    verdict = payload["verdict"]
    assert verdict["passed"] is True
    assert verdict["classification"]["status"] in ("agrees", "tessellation-deficit")
    # the analytic reference is the point of this case, and the engines have to land on it
    reference = payload["case"]["reference"]["volume_m3"]
    got = payload["case"]["routes"]["engine_a_trimesh"]
    assert abs(got - reference) / reference < xv.ACCEPTANCE_TOLERANCE
    assert len(payload["case"]["deficits"]) >= 2


def test_console_script_reports_an_honest_fail_at_a_coarse_resolution(capsys):
    """61 x 13 stations/levels is 0.33 % low: the run must exit 1, not 0."""
    from tensorcad.cli import main as cli_main

    code = cli_main(
        ["volume-xval", "wigley", "--stations", "61", "--levels", "13", "--json"]
    )
    payload = json.loads(capsys.readouterr().out)
    assert code == 1
    gates = {gate["name"]: gate for gate in payload["verdict"]["gates"]}
    assert gates["engines_agree"]["passed"] is True
    assert gates["reference_error"]["passed"] is False
    # and the failure is attributed to tessellation, measured on the ladder
    assert payload["verdict"]["classification"]["status"] == "tessellation-deficit"
    assert 1.5 <= payload["verdict"]["classification"]["observed_order"] <= 2.6


def test_console_script_writes_the_json_it_prints(tmp_path: Path, capsys):
    from tensorcad.cli import main as cli_main

    out = tmp_path / "volume_xval.json"
    code = cli_main(["volume-xval", "wigley", "--stations", "21", "--levels", "9", "--out", str(out)])
    printed = capsys.readouterr().out
    assert out.is_file()
    payload = json.loads(out.read_text())
    assert payload["schema"].startswith("tensorcad.volume_xval_cli/")
    assert "route" in printed and code in (0, 1)


def test_console_script_refuses_an_empty_command():
    from tensorcad.cli import main as cli_main

    assert cli_main(["volume-xval"]) == 2


def test_console_script_reads_an_stl_the_operator_already_has(tmp_path: Path, capsys):
    from tensorcad.cli import main as cli_main

    path = _box_stl(tmp_path / "box.stl")
    code = cli_main(["volume-xval", "mesh", str(path), "--json", "--reference-volume", "30.0"])
    payload = json.loads(capsys.readouterr().out)
    assert code == 0
    routes = payload["case"]["routes"]
    assert routes["engine_a_trimesh"] == pytest.approx(30.0, rel=1e-6)
    assert routes["engine_b_numpy_stl"] == pytest.approx(30.0, rel=1e-6)
    assert payload["verdict"]["classification"]["status"] == "agrees"


def test_console_script_refuses_a_missing_file(tmp_path: Path, capsys):
    from tensorcad.cli import main as cli_main

    code = cli_main(["volume-xval", "mesh", str(tmp_path / "nope.stl")])
    assert code == 2
    assert "refused" in capsys.readouterr().err


# ---------------------------------------------------------------------------
# integration caller: the fleet gate
# ---------------------------------------------------------------------------
def test_fleet_gate_calls_the_independent_engines(tmp_path: Path):
    solid, analytic = _wigley_solid()
    _write_binary_stl(tmp_path / "Wigley_Hull.stl", solid)
    report = audit_fleet(tmp_path, published_volumes={"Wigley": analytic})
    grades = {ship["ship"]: ship for ship in report["ships"]}
    assert grades["Wigley"]["grade"] == VERIFIED
    evidence = report["cross_engine"]["Wigley"]
    assert evidence["status"] == "agrees"
    assert set(evidence["routes"]) == {"divergence", "trimesh", "numpy_stl"}
    assert evidence["engine_spread"] <= xv.ENGINE_AGREEMENT_TOLERANCE
    # the gate's own number and the third-party engine agree, which is the check
    assert abs(evidence["gate_volume_vs_engine_a"]) <= 10 * xv.ENGINE_AGREEMENT_TOLERANCE
    assert evidence["withheld"] == []


def test_fleet_gate_can_switch_the_engines_off(tmp_path: Path):
    solid, analytic = _wigley_solid()
    _write_binary_stl(tmp_path / "Wigley_Hull.stl", solid)
    report = audit_fleet(
        tmp_path, published_volumes={"Wigley": analytic}, cross_engine=False
    )
    assert report["cross_engine"] == {}


def test_fleet_gate_records_a_missing_engine_instead_of_passing_silently(
    tmp_path: Path, monkeypatch
):
    solid, analytic = _wigley_solid(stations=41, levels=9)
    _write_binary_stl(tmp_path / "Wigley_Hull.stl", solid)

    def unavailable(_mesh):
        raise xv.VolumeXvalError("engine B (numpy-stl) is not installed in this environment")

    monkeypatch.setattr(xv, "volume_numpy_stl", unavailable)
    report = audit_fleet(tmp_path, published_volumes={"Wigley": analytic})
    evidence = report["cross_engine"]["Wigley"]
    assert evidence["status"] == "unavailable"
    assert any("numpy-stl" in reason for reason in evidence["withheld"])
    assert "not cross-checked" in evidence["note"]


def test_fleet_gate_downgrades_a_volume_the_engines_cannot_reproduce(
    tmp_path: Path, monkeypatch
):
    """Two engines reading the same triangles must agree; if they do not, refuse."""
    solid, analytic = _wigley_solid(stations=81, levels=17)
    _write_binary_stl(tmp_path / "Wigley_Hull.stl", solid)

    monkeypatch.setattr(xv, "volume_numpy_stl", lambda _mesh: 0.9 * xv.volume_trimesh(_mesh))
    report = audit_fleet(tmp_path, published_volumes={"Wigley": analytic})
    grades = {ship["ship"]: ship for ship in report["ships"]}
    assert grades["Wigley"]["grade"] != VERIFIED
    assert any("not reproducible" in reason for reason in grades["Wigley"]["reasons"])
    assert report["cross_engine"]["Wigley"]["status"] == "engine-disagreement"