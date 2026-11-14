"""The Series 60 family through the product surface: ``tensorcad-hull series60``.

The closed forms are covered elsewhere (``test_series60*.py``).  What is covered
here is the *entry*: that the command a user actually runs reaches the family,
that the exported document carries the family's own analytic reference so a
consumer can recompute every error the entry claims, that the gate refuses a
grid that misses that reference instead of writing a file anyway, and that the
published comparison labels what it proves (closure vs prediction) instead of
calling everything a verification.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

import tensorcad
from tensorcad.cli import main
from tensorcad.series60_cli import main as series60_main


def _run(tmp_path: Path, *extra: str, out: str | None = None):
    output = out or str(tmp_path / "hull")
    argv = ["parent", "--length", "10", "--beam", "1.5", "--draft", "0.6",
            "--cb", "0.70", "--cwp", "0.753", "--lcb", "0.490",
            "--output", output, *extra]
    return main(["series60", *argv])


def test_package_exports_the_family():
    """The kernel has to be reachable from the package, not only from a submodule."""
    assert tensorcad.PARENT_CB == (0.60, 0.65, 0.70, 0.75, 0.80)
    model = tensorcad.generate_series60(tensorcad.Series60Parameters(cb=0.70, cwp=0.753, lcb=0.490))
    assert model.parameters.cb == pytest.approx(0.70)


def test_cli_generates_the_family_and_carries_its_reference(tmp_path):
    code = _run(tmp_path, "--gate", "--stl")
    assert code == 0
    document = json.loads((tmp_path / "hull" / "metadata.json").read_text())

    assert document["model_type"] == "series_60_parametric_displacement_hull"
    assert document["parameters"]["cb"] == pytest.approx(0.70)
    # The reference travels with the artifact: a consumer that reads only the file
    # can recompute every error below without importing this package.
    reference = document["analytic_reference"]
    assert reference["displacement_volume"] == pytest.approx(10.0 * 1.5 * 0.6 * 0.70)

    checked = 0
    for key, value in document["analytic_relative_errors"].items():
        if value is None:  # no closed form for the wetted surface on this family
            continue
        assert abs(value) < 0.02, (key, value)
        checked += 1
    assert checked >= 8

    assert document["gate"]["passed"] is True
    assert document["discretisation"]["within_tolerance"] is True
    for name in ("hull.obj", "offsets.csv", "hull.stl", "metadata.json"):
        path = tmp_path / "hull" / name
        assert path.is_file() and path.stat().st_size > 0


def test_cli_honours_the_block_coefficient(tmp_path):
    """A CLI that silently drops ``--cb`` would still produce a valid-looking hull."""
    for cb in (0.60, 0.75):
        out = tmp_path / f"cb{cb}"
        # No Cwp/LCB here: those are the fidelity targets of one parent form, and a
        # 0.753 waterplane on a Cb = 0.60 hull is an invalid combination by design.
        code = main(["series60", "parent", "--cb", f"{cb}", "--output", str(out)])
        assert code == 0
        document = json.loads((out / "metadata.json").read_text())
        assert document["parameters"]["cb"] == pytest.approx(cb)
        assert document["analytic_reference"]["block_coefficient"] == pytest.approx(cb, abs=1e-12)
        measured = document["hydrostatics"]["block_coefficient"]
        assert abs(measured - cb) / cb < 0.02


def test_cli_refuses_a_grid_that_misses_the_closed_form(tmp_path):
    """The gate has teeth: three vertical levels cannot integrate a zeta**(1/m) section."""
    out = tmp_path / "coarse"
    code = _run(tmp_path, "--nz", "3", "--gate", out=str(out))
    assert code == 1
    assert not (out / "hull.obj").exists()
    assert not (out / "metadata.json").exists()


def test_cli_rejects_invalid_combinations(tmp_path):
    # Cwp without the volume to support it: the deck plan is too full for Cb.
    code = main(["series60", "parent", "--cb", "0.60", "--cwp", "0.95", "--json",
                 "--output", str(tmp_path / "bad")])
    assert code == 2
    assert not (tmp_path / "bad").exists()
    # A bulb needs both halves of its parameterisation.
    code = main(["series60", "parent", "--bulb-length-ratio", "0.02",
                 "--output", str(tmp_path / "bad2")])
    assert code == 2
    # A block coefficient outside the family's calibrated band.
    code = main(["series60", "parent", "--cb", "0.42",
                 "--output", str(tmp_path / "bad3")])
    assert code == 2


def test_family_entry_compares_the_five_parents(tmp_path, monkeypatch):
    """Five parents, published table read not retyped, roles stated honestly."""
    repo = Path(__file__).resolve().parents[1]
    published = repo / "docs" / "data" / "series60_published_parents.json"
    code = series60_main(["family", "--nx", "81", "--nz", "21", "--json",
                          "--published", str(published),
                          "--out", str(tmp_path / "family")])
    assert code == 0
    document = json.loads((tmp_path / "family" / "five_cb_gradient.json").read_text())
    assert [row["cb_requested"] for row in document["rows"]] == [0.60, 0.65, 0.70, 0.75, 0.80]
    assert document["summary"]["within_published_tolerance"] is True
    assert abs(document["summary"]["worst_published_relative"]) < 0.03

    # The role labels are the honest part: Cwp and LCB are targets handed to the
    # generator, so their agreement is a closure check, and Cp is the published
    # ladder the exponents are read from -- calling any of them "independent"
    # would turn a closure into a verification claim.
    roles = {}
    for row in document["rows"]:
        for source, entries in row["error_vs_published"].items():
            for name, entry in entries.items():
                roles.setdefault(name, set()).add(entry["role"])
    assert roles["waterplane_coefficient"] == {"input target (closure check)"}
    assert roles["lcb_over_length"] == {"input target (closure check)"}
    assert roles["prismatic_coefficient"] == {"published ladder (input, not a prediction)"}
    assert roles["block_coefficient"] == {"derived"}


def test_family_entry_without_a_published_table_is_closure_only(tmp_path):
    """No table, no comparison: the entry may not invent one."""
    code = series60_main(["family", "--cb", "0.70", "--nx", "41", "--nz", "13", "--json",
                          "--out", str(tmp_path / "solo")])
    assert code == 0
    document = json.loads((tmp_path / "solo" / "five_cb_gradient.json").read_text())
    assert len(document["rows"]) == 1
    assert document["rows"][0]["error_vs_published"] == {}
    assert document["summary"]["within_published_tolerance"] is True