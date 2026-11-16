"""The Gmsh source differential must be able to FAIL.

The differential's value is entirely in its gate: a run that compares the library with
Gmsh's own arithmetic and then cannot report a disagreement is worse than no run, because
it reads as evidence. These tests pin the three properties that make the arm honest
without needing a Gmsh source tree or a compiled oracle:

* the gate arithmetic counts disagreements (a wholly disagreeing oracle exits 1);
* an absent oracle is exit 2 under ``--require-oracle``, never 0;
* the *build* defines ``HAVE_MESH``. Every Gmsh shape measure is
  ``#if defined(HAVE_MESH) ... #else return 0.;``, so a config header that leaves it
  undefined yields an oracle that answers 0.0 for everything and exits 0 -- the exact
  defect this file guards against.

The compiled-oracle comparison itself is not run here: it needs a Gmsh 4.15.2 source
tree and ~90 s of compiling, so it lives in ``docs/mesh_quality_gmsh_have_mesh_switch.md``
with the numbers it produced.
"""
from __future__ import annotations

import json
import stat
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
DRIVER = REPO / "tools" / "mesh_quality_gmsh_source_diff.py"
BUILD = REPO / "tools" / "mesh_quality_oracles" / "build_gmsh_source_oracle.sh"
HARNESS = REPO / "tools" / "mesh_quality_oracles" / "gmsh_source_oracle.cpp"


def _fake_oracle(path: Path, scale: float) -> Path:
    """An oracle that speaks the ``W`` protocol and returns ``scale`` for every measure.

    ``scale=0.0`` is a deliberate lie: the corpora contain elements whose true measures are
    far from 0, so a driver that cannot see the lie is not gating anything.
    """
    path.write_text(
        "#!/usr/bin/env python3\n"
        "import sys\n"
        "SCALE = " + repr(float(scale)) + "\n"
        "KINDS = {'tri': ('gamma', 'eta', 'minSICN'),\n"
        "         'quad': ('eta', 'area', 'minSICN'),\n"
        "         'tet': ('gamma', 'eta', 'volume'), 'hex': ('minSICN',)}\n"
        "f = open(sys.argv[1])\n"
        "for line in f:\n"
        "    t = line.split()\n"
        "    if not t or t[0] == '#':\n"
        "        continue\n"
        "    if t[0] == 'nodes':\n"
        "        for _ in range(int(t[1])):\n"
        "            next(f)\n"
        "        continue\n"
        "    n = int(t[1])\n"
        "    for _ in range(n):\n"
        "        next(f)\n"
        "    for m in KINDS.get(t[0], ()):\n"
        "        for i in range(n):\n"
        "            print('W %s %s %d %.17g' % (t[0], m, i, SCALE))\n"
    )
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    return path


def test_gate_fails_when_gmsh_arm_disagrees(tmp_path):
    """A wholly disagreeing oracle must exit 1 and be counted, element by element."""
    oracle = _fake_oracle(tmp_path / "zero_oracle", 0.0)
    artifact = tmp_path / "diff.json"
    proc = subprocess.run(
        [sys.executable, str(DRIVER), "--oracle", str(oracle),
         "--require-oracle", "--artifact", str(artifact)],
        capture_output=True, text=True,
    )
    assert proc.returncode == 1, proc.stdout + proc.stderr
    report = json.loads(artifact.read_text())
    assert report["pairs_compared"] >= 8          # tri, quad, tet over five meshes
    assert report["elements_beyond_tolerance"] > 0
    assert "pairs beyond" not in proc.stdout      # stale wording would mean stale code
    # every compared metric is reported as DIFFER, none as the old phantom "COMPARED"
    statuses = {info["status"] for entry in report["meshes"].values()
                for info in entry.get("metrics", {}).values()}
    assert statuses == {"DIFFER"}, statuses


def test_absent_oracle_is_exit_2_not_success(tmp_path):
    """No oracle must never look like agreement."""
    artifact = tmp_path / "diff.json"
    proc = subprocess.run(
        [sys.executable, str(DRIVER), "--oracle", str(tmp_path / "does_not_exist"),
         "--require-oracle", "--artifact", str(artifact)],
        capture_output=True, text=True,
    )
    assert proc.returncode == 2, proc.stdout + proc.stderr
    report = json.loads(artifact.read_text())
    assert report["oracle_ran"] is False
    assert report["pairs_compared"] == 0


def test_build_defines_have_mesh():
    """`HAVE_MESH` is the switch every Gmsh shape measure is guarded by.

    Undefined -> ``#else return 0.;`` -> an oracle that reports 0.0 for gamma/eta/minSICN
    and exits 0. That is the defect this build script now cannot reintroduce silently.
    """
    script = BUILD.read_text()
    assert 'DEFINE_ON="HAVE_MESH"' in script
    assert "open(sys.argv[1]).read()" in script, "config generator must read the template"
    assert "'#define %s 1' % name if name in define_on" in script


def test_harness_emits_the_measured_metrics():
    """The harness must ask for the measures the driver maps, and for hex at all."""
    cpp = HARNESS.read_text()
    for token in ('emit("tri", "gamma"', 'emit("tet", "volume"',
                  'emit("tet", "minSICN"', 'MHexahedron el(',
                  'emit("hex", "minSIGE"'):
        assert token in cpp, token


def test_oracle_protocol_lines_are_parsed_by_the_shared_reader(tmp_path):
    """A fake oracle's lines round-trip through the shared parser (guards the protocol)."""
    sys.path.insert(0, str(REPO / "tools"))
    import mesh_quality_external_diff as drv

    mesh = tmp_path / "m.mkm"
    drv.write_mesh_file(mesh, "tri", [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)],
                        [[0, 1, 2]])
    oracle = _fake_oracle(tmp_path / "zero_oracle", 0.0)
    status, data = drv.run_oracle(str(oracle), mesh, "W")
    assert status == "RAN"
    assert data["tri"]["gamma"] == [0.0]


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(pytest.main([__file__, "-q"]))